import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from scipy.optimize import brentq


def load_data(metocean_path, fls_path):
    """Load FLS proxy results and merge with metocean data."""
    met_df = pd.read_excel(metocean_path, sheet_name='metocean')
    fls_df = pd.read_excel(fls_path)
    # Tp already in FLS output; only merge Hs and wave_dir from metocean
    merge_cols = ['time', 'Hs', 'wave_dir']
    if 'Tp' not in fls_df.columns:
        merge_cols.append('Tp')
    df = pd.merge(fls_df, met_df[merge_cols], on='time')
    return df


def plot_frequency(df, hs_bins=10, tp_bins=10, graphs_dir='graphs'):
    """Plot Hs-Tp bin frequency."""
    hs_edges = np.linspace(df['Hs'].min(), df['Hs'].max(), hs_bins + 1)
    tp_edges = np.linspace(df['Tp'].min(), df['Tp'].max(), tp_bins + 1)

    H = np.zeros((tp_bins, hs_bins))
    for i in range(hs_bins):
        for j in range(tp_bins):
            mask = (
                (df['Hs'] >= hs_edges[i]) & (df['Hs'] < hs_edges[i+1]) &
                (df['Tp'] >= tp_edges[j]) & (df['Tp'] < tp_edges[j+1])
            )
            H[j, i] = mask.sum()

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(H, origin='lower', aspect='auto',
                   extent=[hs_edges[0], hs_edges[-1], tp_edges[0], tp_edges[-1]],
                   norm=PowerNorm(gamma=0.4), cmap='plasma')
    ax.set_xlabel('Hs [m]')
    ax.set_ylabel('Tp [s]')
    ax.set_title('Bin Frequency (count per Hs-Tp bin)')
    plt.colorbar(im, ax=ax, label='Count')
    plt.tight_layout()
    plt.savefig(os.path.join(graphs_dir, 'pendulum_frequency.png'), dpi=150)
    plt.close()


def plot_heatmaps(df, hs_bins=10, tp_bins=10, graphs_dir='graphs'):
    """Plot Hs-Tp heatmaps of key pendulum model outputs."""
    variables = [
        ('Moment_Crest_kNm',  'Moment Crest [kNm]'),
        ('Moment_Trough_kNm', 'Moment Trough [kNm]'),
        ('Strain_Range',      'Strain Range [-]'),
        ('FLS_Proxy',         'FLS Proxy [-]'),
    ]

    hs_edges = np.linspace(df['Hs'].min(), df['Hs'].max(), hs_bins + 1)
    tp_edges = np.linspace(df['Tp'].min(), df['Tp'].max(), tp_bins + 1)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for ax, (var, label) in zip(axes, variables):
        H = np.zeros((tp_bins, hs_bins))
        for i in range(hs_bins):
            for j in range(tp_bins):
                mask = (
                    (df['Hs'] >= hs_edges[i]) & (df['Hs'] < hs_edges[i+1]) &
                    (df['Tp'] >= tp_edges[j]) & (df['Tp'] < tp_edges[j+1])
                )
                H[j, i] = df.loc[mask, var].sum()

        im = ax.imshow(H, origin='lower', aspect='auto',
                       extent=[hs_edges[0], hs_edges[-1], tp_edges[0], tp_edges[-1]],
                       norm=PowerNorm(gamma=0.4), cmap='plasma')
        ax.set_xlabel('Hs [m]')
        ax.set_ylabel('Tp [s]')
        ax.set_title(f'Cumulative {label}')
        plt.colorbar(im, ax=ax)

    plt.tight_layout()
    plt.savefig(os.path.join(graphs_dir, 'pendulum_heatmaps.png'), dpi=150)
    plt.close()


def plot_damage_roses(df, n_sectors=12, graphs_dir='graphs'):
    """Plot directional damage roses for key parameters."""
    variables = [
        ('Moment_Crest_kNm',  'Moment Crest [kNm]'),
        ('Moment_Trough_kNm', 'Moment Trough [kNm]'),
        ('Strain_Range',      'Strain Range [-]'),
        ('FLS_Proxy',         'FLS Proxy [-]'),
    ]

    sector_edges = np.linspace(0, 360, n_sectors + 1)
    sector_centers = (sector_edges[:-1] + sector_edges[1:]) / 2
    theta = np.radians(sector_centers)
    width = np.radians(360 / n_sectors)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), subplot_kw={'projection': 'polar'})
    axes = axes.flatten()

    for ax, (var, label) in zip(axes, variables):
        radii = []
        for i in range(n_sectors):
            mask = (df['wave_dir'] >= sector_edges[i]) & (df['wave_dir'] < sector_edges[i+1])
            radii.append(df.loc[mask, var].sum())

        ax.bar(theta, radii, width=width, alpha=0.7)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_title(f'{label}')

    plt.tight_layout()
    plt.savefig(os.path.join(graphs_dir, 'pendulum_damage_roses.png'), dpi=150)
    plt.close()


def _N_from_epsilon_scalar(eps, A_el, b_el, A_pl, b_pl):
    """Solve epsilon = A_el * N^b_el + A_pl * N^b_pl for N."""
    def residual(logN):
        N = 10**logN
        return A_el * N**b_el + A_pl * N**b_pl - eps

    try:
        return 10**brentq(residual, -1, 15)
    except ValueError:
        return np.nan


def plot_strain_n_curve(df, graphs_dir='graphs'):
    """
    Plot epsilon-N fatigue curves (P2169 & P2159) with simulated
    strain values from the pendulum model overlaid.
    """
    # ── Reference curve parameters ──
    # P2169 (Design Report 10.3.4 Lead sheath)
    A_el, b_el = 0.0025, -0.18
    A_pl, b_pl = 66.0, -1.153

    # P2159 (client-provided)
    logA, m_sn = -6.2369, -3.8849

    N_ref = np.logspace(0, 9, 200)

    # Curve 1: P2169 direct  epsilon(N) = A_el*N^b_el + A_pl*N^b_pl
    eps_total = A_el * N_ref**b_el + A_pl * N_ref**b_pl

    # Curve 2: P2169 validation (inverse)
    eps_range = np.logspace(-5, 1, 100)
    N_validation = np.array([
        _N_from_epsilon_scalar(e, A_el, b_el, A_pl, b_pl) for e in eps_range
    ])

    # Curve 3: P2159  epsilon(N) = 10^((log10(N) - logA) / m)
    eps_p2159 = 10**((np.log10(N_ref) - logA) / m_sn)

    # ── Simulated data (cumulative cycles per strain level) ──
    df_sim = df.loc[(df['Strain_Range'] > 1e-8) & (df['Tp'] > 0)].copy()
    df_sim['n_cycles'] = 3600.0 / df_sim['Tp']

    # ── Plot ──
    fig, ax = plt.subplots(figsize=(10, 7))

    ax.loglog(N_ref, eps_total, 'k-', linewidth=2,
              label=r'$\epsilon_a$-N P2169')
    ax.loglog(N_validation, eps_range, 'r--', linewidth=2,
              label=r'$\epsilon_a$-N P2169 Validation')
    ax.loglog(N_ref, eps_p2159, 'b-.', linewidth=1.5,
              label='P2159')

    if len(df_sim) > 0:
        # Bin by strain range (log-spaced) and sum total cycles per bin
        log_min = np.log10(df_sim['Strain_Range'].min())
        log_max = np.log10(df_sim['Strain_Range'].max())
        bins = np.logspace(log_min, log_max, 201)
        df_sim['strain_bin'] = pd.cut(df_sim['Strain_Range'], bins=bins)

        agg_dict = {'n_cycles': 'sum', 'Strain_Range': 'mean'}
        if 'Hs' in df_sim.columns:
            agg_dict['Hs'] = 'mean'

        grouped = df_sim.groupby('strain_bin', observed=True).agg(agg_dict).dropna()

        strain_agg = grouped['Strain_Range'].values
        cycles_agg = grouped['n_cycles'].values

        if 'Hs' in grouped.columns:
            sc = ax.scatter(cycles_agg, strain_agg,
                            c=grouped['Hs'].values, s=12, alpha=0.6,
                            cmap='coolwarm', zorder=5,
                            label='Simulated (cumulative)')
            plt.colorbar(sc, ax=ax, label='Hs [m]', pad=0.02)
        else:
            ax.scatter(cycles_agg, strain_agg,
                       c='green', s=12, alpha=0.6, zorder=5,
                       label='Simulated (cumulative)')

    ax.set_xlabel('Number of cycles, N', fontsize=12)
    ax.set_ylabel(r'Strain amplitude, $\epsilon_a$ [mm/mm]', fontsize=12)
    ax.set_title(r'$\epsilon$-N Fatigue Curve with Pendulum Model Results', fontsize=14)
    ax.set_xlim(1, 1e9)
    ax.set_ylim(1e-8, 1e1)
    ax.legend(loc='upper right')
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(graphs_dir, 'strain_n_curve.png'), dpi=150)
    plt.close()


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    graphs_dir = os.path.join(script_dir, "graphs")
    os.makedirs(graphs_dir, exist_ok=True)

    metocean_file = os.path.join(script_dir, "metocean.xlsx")
    fls_file = os.path.join(script_dir, "fls_proxy_range_method.xlsx")

    hs_bins = 50
    tp_bins = 50
    n_sectors = 36

    df = load_data(metocean_file, fls_file)
    print(f"Loaded {len(df)} rows")
    print(df.head())

    plot_frequency(df, hs_bins, tp_bins, graphs_dir)
    print("Saved: pendulum_frequency.png")

    plot_heatmaps(df, hs_bins, tp_bins, graphs_dir)
    print("Saved: pendulum_heatmaps.png")

    plot_damage_roses(df, n_sectors, graphs_dir)
    print("Saved: pendulum_damage_roses.png")

    plot_strain_n_curve(df, graphs_dir)
    print("Saved: strain_n_curve.png")
