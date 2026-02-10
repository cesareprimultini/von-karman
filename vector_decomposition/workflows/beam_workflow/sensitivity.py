# global sensitivity analysis using Sobol indices on the actual beam model.
# answers: "which metocean variable drives FLS the most, and by how much?"
#
# Sobol decomposes the total FLS variance into fractions per input:
#   S1 (first-order) = direct effect of that variable alone
#   ST (total-order) = direct effect + all interactions with other variables
#   ST - S1          = importance through interactions only
#
# all evaluations go through the real physics (Airy + beam), no surrogate.
# requires: pip install SALib

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from SALib.sample import saltelli
from SALib.analyze import sobol
from matplotlib.colors import LogNorm

from airy_wave_velocities import compute_wave_velocities
from fls_proxy_ff_v2 import (build_segments_3d, _solve_beam_force_range,
                              get_curvature_from_moment)

# --- constants (match fls_proxy_ff_v2.py) ---
CPS_HEADING = 65
D_PIPE = 0.251
RHO = 1025
CD = 1.2
WATER_DEPTH = 53.5
M_FATIGUE = 3

PARAM_NAMES = ['Hs', 'Tp', 'wave_dir', 'current_dir', 'Uc', 'water_level']


def _setup_beam(segments_path):
    """Load segment geometry and precompute all fixed beam parameters."""
    seg_df = pd.read_excel(segments_path)
    segments, bellmouth, touchdown = build_segments_3d(seg_df, CPS_HEADING)

    beam_vec = touchdown - bellmouth
    L_beam = np.linalg.norm(beam_vec)
    beam_axis = beam_vec / L_beam

    L_arc = sum(seg['length'] for seg in segments)

    # Arc-length based load positions (matches fls_proxy_ff_v2.py)
    load_positions = []
    arc_running = 0.0
    for seg in segments:
        arc_to_mid = arc_running + seg['length'] / 2.0
        a_i = (arc_to_mid / L_arc) * L_beam
        load_positions.append(max(0.0, min(a_i, L_beam)))
        arc_running += seg['length']

    # Arc-length based edge positions (segment boundaries)
    edge_positions = [0.0]
    arc_running = 0.0
    for seg in segments:
        arc_running += seg['length']
        ep = (arc_running / L_arc) * L_beam
        edge_positions.append(max(0.0, min(ep, L_beam)))

    midspan = L_beam / 2
    eval_points = sorted(set(edge_positions + load_positions + [midspan]))

    return dict(
        segments=segments, L_beam=L_beam, beam_axis=beam_axis,
        load_positions=load_positions, edge_positions=edge_positions,
        eval_points=eval_points, midspan=midspan,
        z_mids=seg_df['z_mid'].values, seg_ids=seg_df['seg_id'].values
    )


def _evaluate_fls(Hs, Tp, wave_dir, current_dir, Uc, water_level, beam):
    """Run Airy + beam model for one parameter set, return FLS proxy."""
    d = WATER_DEPTH + water_level

    # u_w per segment via Airy wave theory
    u_w_all = {}
    for z_mid, seg_id in zip(beam['z_mids'], beam['seg_ids']):
        u_w, _, _, _ = compute_wave_velocities(Hs, Tp, d, z_mid)
        u_w_all[seg_id] = np.array([u_w])

    wave_rad = np.radians(wave_dir)
    wave_dir_vec = np.array([np.sin(wave_rad), np.cos(wave_rad), 0])

    c_rad = np.radians(current_dir)
    curr_vec = Uc * np.array([np.sin(c_rad), np.cos(c_rad), 0])

    rng = _solve_beam_force_range(
        segments=beam['segments'], u_w_all=u_w_all, idx=0,
        wave_dir_vec=wave_dir_vec, curr_vec=curr_vec,
        beam_axis=beam['beam_axis'],
        load_positions=beam['load_positions'],
        edge_positions=beam['edge_positions'],
        eval_points=beam['eval_points'],
        midspan=beam['midspan'], L_beam=beam['L_beam'],
        rho=RHO, Cd=CD, D=D_PIPE
    )

    strain_range = get_curvature_from_moment(rng['M_range_max']) * (D_PIPE / 2.0)
    n_cycles = 3600.0 / Tp
    return n_cycles * strain_range ** M_FATIGUE


def _get_bounds(metocean_path):
    """Parameter ranges from the dataset (1st-99th percentile)."""
    met = pd.read_excel(metocean_path, sheet_name='metocean')
    bounds = []
    for col in PARAM_NAMES:
        lo = met[col].quantile(0) # min
        hi = met[col].quantile(1) # max (can be a quantile, e.g. 0.99)
        bounds.append([lo, hi])

    # directions: full circle so Sobol can explore all orientations
    bounds[2] = [0.0, 360.0]   # wave_dir
    bounds[3] = [0.0, 360.0]   # current_dir

    # safety floors
    bounds[0][0] = max(bounds[0][0], 0.05)  # Hs > 0
    bounds[1][0] = max(bounds[1][0], 1.0)   # Tp > 1s
    bounds[4][0] = max(bounds[4][0], 0.0)   # Uc >= 0

    return bounds


def run_sobol(segments_path, metocean_path, N=1024):
    """Generate Saltelli samples, evaluate the beam model, compute Sobol indices."""
    beam = _setup_beam(segments_path)
    bounds = _get_bounds(metocean_path)

    problem = {
        'num_vars': len(PARAM_NAMES),
        'names': PARAM_NAMES,
        'bounds': bounds
    }

    samples = saltelli.sample(problem, N)
    n_samples = len(samples)
    print(f"Sobol analysis: {n_samples} samples ({N} base, {len(PARAM_NAMES)} params)")
    print(f"Bounds:")
    for name, (lo, hi) in zip(PARAM_NAMES, bounds):
        print(f"  {name:<15} [{lo:.3f}, {hi:.3f}]")

    Y = np.zeros(n_samples)
    for i in range(n_samples):
        Y[i] = _evaluate_fls(*samples[i], beam)
        if (i + 1) % 2000 == 0:
            print(f"  {i + 1} / {n_samples}")

    Si = sobol.analyze(problem, Y)
    return Si


def _get_medians(metocean_path):
    """Median value for each parameter (used as baseline for 2D sweeps)."""
    met = pd.read_excel(metocean_path, sheet_name='metocean')
    return {col: met[col].median() for col in PARAM_NAMES}


def plot_sobol(Si, graphs_dir):
    """Bar chart of first-order and total-order Sobol indices."""
    S1 = Si['S1']
    ST = Si['ST']

    order = np.argsort(ST)[::-1]
    names_sorted = [PARAM_NAMES[i] for i in order]

    x = np.arange(len(PARAM_NAMES))
    w = 0.35

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(x + w / 2, S1[order], w, xerr=Si['S1_conf'][order],
            label='First-order (S1)', color='steelblue',
            edgecolor='k', linewidth=0.4)
    ax.barh(x - w / 2, ST[order], w, xerr=Si['ST_conf'][order],
            label='Total-order (ST)', color='coral',
            edgecolor='k', linewidth=0.4)

    ax.set_yticks(x)
    ax.set_yticklabels(names_sorted)
    ax.set_xlabel('Sobol index (fraction of FLS variance)')
    ax.set_title('Global Sensitivity - Sobol Indices')
    ax.legend(loc='lower right')
    ax.set_xlim(left=-0.05)

    plt.tight_layout()
    out = os.path.join(graphs_dir, 'sobol_sensitivity.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}")


def plot_heatmaps(beam, medians, bounds, graphs_dir, n_grid=80):
    """2D heatmaps of FLS for key parameter pairs, all others held at median."""

    # pairs that capture the most important interactions
    pairs = [
        ('Hs',          'wave_dir',     'Wave height vs direction'),
        ('Uc',          'current_dir',  'Current magnitude vs direction'),
        ('wave_dir',    'current_dir',  'Wave vs current direction'),
        ('Hs',          'Uc',           'Wave height vs current magnitude'),
        ('Hs',          'Tp',           'Wave height vs period'),
        ('Uc',          'wave_dir',     'Current magnitude vs wave direction'),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    for ax, (p1, p2, title) in zip(axes.flatten(), pairs):
        i1 = PARAM_NAMES.index(p1)
        i2 = PARAM_NAMES.index(p2)

        v1 = np.linspace(bounds[i1][0], bounds[i1][1], n_grid)
        v2 = np.linspace(bounds[i2][0], bounds[i2][1], n_grid)
        Z = np.zeros((n_grid, n_grid))

        for j, val2 in enumerate(v2):
            for k, val1 in enumerate(v1):
                params = dict(medians)
                params[p1] = val1
                params[p2] = val2
                Z[j, k] = _evaluate_fls(
                    params['Hs'], params['Tp'], params['wave_dir'],
                    params['current_dir'], params['Uc'], params['water_level'],
                    beam
                )

        im = ax.pcolormesh(v1, v2, Z, shading='auto', cmap='viridis', norm=LogNorm(vmin=max(Z[Z > 0].min(), 1e-30)))
        fig.colorbar(im, ax=ax, label='FLS proxy')
        ax.set_xlabel(p1)
        ax.set_ylabel(p2)
        ax.set_title(title)

        # mark the pipe heading on direction axes
        for param, axis_fn in [(p1, ax.axvline), (p2, ax.axhline)]:
            if 'dir' in param:
                axis_fn(CPS_HEADING, color='cyan', ls='--', lw=1, label=f'pipe heading ({CPS_HEADING} deg)')
                axis_fn((CPS_HEADING + 180) % 360, color='cyan', ls='--', lw=1)
                # perpendicular to pipe
                axis_fn((CPS_HEADING + 90) % 360, color='lime', ls=':', lw=1, label=f'perp to pipe')
                axis_fn((CPS_HEADING + 270) % 360, color='lime', ls=':', lw=1)

# FLS sensitivity heatmaps (other params at median)
# cyan dashed = inline with pipe, green dotted = perpendicular to pipe

    plt.title('FLS sensitivity heatmaps')
    plt.tight_layout()
    out = os.path.join(graphs_dir, 'sensitivity_heatmaps.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}")


def plot_cross_correlation(metocean_path, graphs_dir):
    """Cross-correlation heatmap of raw metocean inputs."""
    met = pd.read_excel(metocean_path, sheet_name='metocean')
    met = met[PARAM_NAMES].dropna()
    corr = met.corr()

    n = len(PARAM_NAMES)
    fig, ax = plt.subplots(figsize=(n * 1.1 + 1, n * 0.9 + 1))
    mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
    sns.heatmap(corr, annot=True, fmt='.2f', cmap='RdBu_r',
                center=0, vmin=-1, vmax=1, linewidths=0.5,
                mask=mask, square=True,
                xticklabels=PARAM_NAMES, yticklabels=PARAM_NAMES, ax=ax)
    ax.set_title('Cross-correlation of metocean inputs')
    plt.tight_layout()
    out = os.path.join(graphs_dir, 'sensitivity_cross_correlation.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}")


def plot_input_distributions(metocean_path, graphs_dir):
    """Violin plots showing distribution of each raw metocean input."""
    met = pd.read_excel(metocean_path, sheet_name='metocean')
    met = met[PARAM_NAMES].dropna()

    fig, axes = plt.subplots(1, len(PARAM_NAMES),
                             figsize=(2.5 * len(PARAM_NAMES), 4), squeeze=False)

    for ax, feat in zip(axes.flatten(), PARAM_NAMES):
        sns.violinplot(y=met[feat], ax=ax, inner='quartile',
                       color=sns.cubehelix_palette(1, rot=-.25, light=.7)[0])
        ax.set_title(feat)
        ax.set_ylabel('')

    plt.suptitle('Distribution of metocean inputs', y=1.02)
    plt.tight_layout()
    out = os.path.join(graphs_dir, 'sensitivity_distributions.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}")


if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))
    graphs_dir = os.path.join(script_dir, 'graphs')
    os.makedirs(graphs_dir, exist_ok=True)

    segments_file = os.path.join(script_dir, 'segments.xlsx')
    metocean_file = os.path.join(script_dir, 'metocean.xlsx')

    # --- Sobol indices ---
    Si = run_sobol(segments_file, metocean_file, N=4096)

    print("\nSobol indices:")
    print(f"{'Parameter':<15} {'S1':>8} {'ST':>8} {'Interaction':>12}")
    print("-" * 45)
    for i, name in enumerate(PARAM_NAMES):
        s1 = Si['S1'][i]
        st = Si['ST'][i]
        print(f"{name:<15} {s1:>8.3f} {st:>8.3f} {st - s1:>12.3f}")

    plot_sobol(Si, graphs_dir)

    # --- 2D heatmaps ---
    beam = _setup_beam(segments_file)
    medians = _get_medians(metocean_file)
    bounds = _get_bounds(metocean_file)

    print(f"\nMedian values used for heatmaps:")
    for name in PARAM_NAMES:
        print(f"  {name:<15} {medians[name]:.3f}")

    print("\nComputing heatmaps...")
    plot_heatmaps(beam, medians, bounds, graphs_dir)

    # --- input data plots ---
    plot_cross_correlation(metocean_file, graphs_dir)
    plot_input_distributions(metocean_file, graphs_dir)
