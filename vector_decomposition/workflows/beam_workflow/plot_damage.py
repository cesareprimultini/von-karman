import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
import seaborn as sns

MOMENT_SCALE = 1e-6   # divide moments by this (show in MN·m)
FLS_SCALE = 1e-12     # divide FLS by this


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
    plt.savefig(os.path.join(graphs_dir, 'frequency.png'), dpi=150)
    plt.close()

def plot_heatmaps(df, hs_bins=10, tp_bins=10, graphs_dir='graphs'):
    variables = [
        ('M_touchdown', f'Moment at TDP [×{MOMENT_SCALE:.0e}]', MOMENT_SCALE),
        ('M_bellmouth', f'Moment at BME [×{MOMENT_SCALE:.0e}]', MOMENT_SCALE),
        ('FLS_Proxy_tdp', f'FLS Proxy at TDP [×{FLS_SCALE:.0e}]', FLS_SCALE),
        ('FLS_Proxy_bme', f'FLS Proxy at BME [×{FLS_SCALE:.0e}]', FLS_SCALE)
    ]

    hs_edges = np.linspace(df['Hs'].min(), df['Hs'].max(), hs_bins + 1)
    tp_edges = np.linspace(df['Tp'].min(), df['Tp'].max(), tp_bins + 1)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for ax, (var, label, scale) in zip(axes, variables):
        H = np.zeros((tp_bins, hs_bins))
        for i in range(hs_bins):
            for j in range(tp_bins):
                mask = (
                    (df['Hs'] >= hs_edges[i]) & (df['Hs'] < hs_edges[i+1]) &
                    (df['Tp'] >= tp_edges[j]) & (df['Tp'] < tp_edges[j+1])
                )
                H[j, i] = df.loc[mask, var].sum() / scale

        im = ax.imshow(H, origin='lower', aspect='auto',
                       extent=[hs_edges[0], hs_edges[-1], tp_edges[0], tp_edges[-1]],
                       norm=PowerNorm(gamma=0.4), cmap='plasma')
        ax.set_xlabel('Hs [m]')
        ax.set_ylabel('Tp [s]')
        ax.set_title(f'Cumulative {label}')
        plt.colorbar(im, ax=ax)

    plt.tight_layout()
    plt.savefig(os.path.join(graphs_dir, 'heatmaps.png'), dpi=150)
    plt.close()

def plot_damage_roses(df, n_sectors=12, graphs_dir='graphs'):
    variables = [
        ('M_touchdown', f'Moment at TDP [×{MOMENT_SCALE:.0e}]', MOMENT_SCALE),
        ('M_bellmouth', f'Moment at BME [×{MOMENT_SCALE:.0e}]', MOMENT_SCALE),
        ('FLS_Proxy_tdp', f'FLS Proxy at TDP [×{FLS_SCALE:.0e}]', FLS_SCALE),
        ('FLS_Proxy_bme', f'FLS Proxy at BME [×{FLS_SCALE:.0e}]', FLS_SCALE)
    ]

    sector_edges = np.linspace(0, 360, n_sectors + 1)
    sector_centers = (sector_edges[:-1] + sector_edges[1:]) / 2
    theta = np.radians(sector_centers)
    width = np.radians(360 / n_sectors)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), subplot_kw={'projection': 'polar'})
    axes = axes.flatten()

    for ax, (var, label, scale) in zip(axes, variables):
        radii = []
        for i in range(n_sectors):
            mask = (df['wave_dir'] >= sector_edges[i]) & (df['wave_dir'] < sector_edges[i+1])
            radii.append(df.loc[mask, var].sum() / scale)

        ax.bar(theta, radii, width=width, alpha=0.7)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_title(f'{label}')

    plt.tight_layout()
    plt.savefig(os.path.join(graphs_dir, 'damage_roses.png'), dpi=150)
    plt.close()


def _rank_bins(df, fls_col, hs_bins, tp_bins):
    """Bin df by Hs-Tp and return (df_with_bin_cols, bin_ranking_df)."""
    hs_edges = np.linspace(df['Hs'].min(), df['Hs'].max(), hs_bins + 1)
    tp_edges = np.linspace(df['Tp'].min(), df['Tp'].max(), tp_bins + 1)
    df = df.copy()
    df['hs_bin'] = pd.cut(df['Hs'], bins=hs_edges, include_lowest=True)
    df['tp_bin'] = pd.cut(df['Tp'], bins=tp_edges, include_lowest=True)
    bin_fls = (df.groupby(['hs_bin', 'tp_bin'], observed=True)[fls_col]
               .sum().reset_index()
               .sort_values(fls_col, ascending=False)
               .reset_index(drop=True))
    return df, bin_fls


def plot_pareto_bins(df, hs_bins=20, tp_bins=20, n_top=30,
                     graphs_dir='graphs'):
    """
    Pareto chart of Hs-Tp bins ranked by cumulative BME FLS proxy,
    with a cumulative-% line and top-N highlighted on the Hs-Tp heatmap.
    """
    fls_col = 'FLS_Proxy_bme' if 'FLS_Proxy_bme' in df.columns else 'FLS_proxy_bme'
    df, bin_fls = _rank_bins(df, fls_col, hs_bins, tp_bins)

    total_fls = bin_fls[fls_col].sum()
    bin_fls['pct'] = 100.0 * bin_fls[fls_col] / total_fls
    bin_fls['cum_pct'] = bin_fls['pct'].cumsum()
    bin_fls['hs_mid'] = bin_fls['hs_bin'].apply(lambda b: b.mid)
    bin_fls['tp_mid'] = bin_fls['tp_bin'].apply(lambda b: b.mid)

    top = bin_fls.head(n_top)
    top_labels = [f"Hs={r.hs_mid:.1f} Tp={r.tp_mid:.1f}"
                  for _, r in top.iterrows()]

    fig, ax1 = plt.subplots(figsize=(max(12, n_top * 0.45), 6))

    colors = sns.cubehelix_palette(n_top, rot=-.25, light=.7)
    ax1.bar(range(n_top), top['pct'], color=colors, edgecolor='k',
            linewidth=0.4)
    ax1.set_xticks(range(n_top))
    ax1.set_xticklabels(top_labels, rotation=70, ha='right', fontsize=7)
    ax1.set_ylabel('Bin contribution to total FLS [%]')
    ax1.set_title(f'Hs-Tp bins ranked by cumulative FLS (BME)')

    ax2 = ax1.twinx()
    ax2.plot(range(n_top), top['cum_pct'].values, 'k-o', markersize=4)
    ax2.set_ylabel('Cumulative [%]')
    ax2.set_ylim(0, 105)
    ax2.axhline(80, color='grey', ls='--', lw=0.8, alpha=0.6)

    plt.tight_layout()
    plt.savefig(os.path.join(graphs_dir, 'ranked_fls_bins.png'),
                dpi=150, bbox_inches='tight')
    plt.close()


def plot_strain_n_curve(df, design_life=25, graphs_dir='graphs'):
    """
    Strain-N curve (P2169 lead sheath) with simulated data aggregated by
    strain level and coloured by Hs.
    """
    # ── Moreno et al. (2021) – strain rate 1E-3 s⁻¹, amplitude × 2 = range ──
    N_ref = np.logspace(0, 9, 300)
    eps_smooth    = 2 * (0.0033 * N_ref**-0.18 + 85.863 * N_ref**-1.153)   # smooth
    eps_irregular = 2 * (0.0025 * N_ref**-0.18 + 66.0   * N_ref**-1.153)   # irregular

    # Viespoli 2020
    eps_viespoli = 9.2593*N_ref**-0.346

    # ── Simulated data: aggregate by strain level ──
    df_sim = df.loc[(df['Strain_Range'] > 1e-8) & (df['Tp'] > 0)].copy()
    hindcast_hours = (pd.to_datetime(df['time']).max() - pd.to_datetime(df['time']).min()).total_seconds() / 3600
    hindcast_years = hindcast_hours / 8766  # hours per year (365.25 days)
    life_scale = design_life / hindcast_years
    df_sim['n_cycles'] = (3600.0 / df_sim['Tp']) * life_scale

    log_min = np.log10(df_sim['Strain_Range'].min())
    log_max = np.log10(df_sim['Strain_Range'].max())
    bins = np.logspace(log_min, log_max, 201)
    df_sim['strain_bin'] = pd.cut(df_sim['Strain_Range'], bins=bins)

    agg_dict = {'n_cycles': 'sum', 'Strain_Range': 'mean'}
    if 'Hs' in df_sim.columns:
        agg_dict['Hs'] = 'mean'

    grouped = (df_sim.groupby('strain_bin', observed=True)
               .agg(agg_dict).dropna())

    strain_agg = grouped['Strain_Range'].values
    cycles_agg = grouped['n_cycles'].values

    # ── Plot ──
    fig, ax = plt.subplots(figsize=(10, 7))

    ax.loglog(N_ref, eps_smooth, 'b-', linewidth=2,
              label='Moreno (2021) Smooth')
    ax.loglog(N_ref, eps_irregular, 'r--', linewidth=2,
              label='Moreno (2021) Irregular')
    ax.loglog(N_ref, eps_viespoli, 'g-.', linewidth=2,
              label='Viespoli (2020)')

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
    ax.set_ylabel(r'Strain range, $\Delta\epsilon$ [mm/mm]', fontsize=12)
    ax.set_title(r'$\epsilon$-N Fatigue Curve',
                 fontsize=13)
    ax.set_xlim(1, 1e9)
    ax.set_ylim(1e-8, 1e1)
    ax.legend(loc='upper right')
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(graphs_dir, 'strain_n_curve.png'),
                dpi=150, bbox_inches='tight')
    plt.close()


def plot_violin_top_bins(df, metocean_path, hs_bins=50, tp_bins=50,
                         n_top=10, graphs_dir='graphs'):
    """
    Half-violin plots with jittered points, median and 95% CI for
    wave_dir, current_dir, water_level, and Uc across the top-N worst
    Hs-Tp bins ranked by cumulative BME FLS proxy.  One figure per variable.
    """
    met_df = pd.read_excel(metocean_path, sheet_name='metocean')
    plot_vars = ['wave_dir', 'current_dir', 'water_level', 'Uc']
    merge_cols = ['time'] + [c for c in plot_vars if c not in df.columns]
    if len(merge_cols) > 1:
        df = pd.merge(df, met_df[merge_cols], on='time', how='left')

    fls_col = 'FLS_Proxy_bme' if 'FLS_Proxy_bme' in df.columns else 'FLS_proxy_bme'

    hs_edges = np.linspace(df['Hs'].min(), df['Hs'].max(), hs_bins + 1)
    tp_edges = np.linspace(df['Tp'].min(), df['Tp'].max(), tp_bins + 1)

    df = df.copy()
    df['hs_bin'] = pd.cut(df['Hs'], bins=hs_edges, include_lowest=True)
    df['tp_bin'] = pd.cut(df['Tp'], bins=tp_edges, include_lowest=True)

    bin_fls = (df.groupby(['hs_bin', 'tp_bin'], observed=True)[fls_col]
               .sum().reset_index()
               .sort_values(fls_col, ascending=False))
    top_bins = bin_fls.head(n_top)

    # Collect data per bin
    bin_labels, bin_data = [], {v: [] for v in plot_vars}
    for rank, (_, row) in enumerate(top_bins.iterrows(), start=1):
        mask = (df['hs_bin'] == row['hs_bin']) & (df['tp_bin'] == row['tp_bin'])
        sub = df.loc[mask]
        hs_mid, tp_mid = row['hs_bin'].mid, row['tp_bin'].mid
        bin_labels.append(f"#{rank}  Hs={hs_mid:.1f} Tp={tp_mid:.1f}")
        for v in plot_vars:
            bin_data[v].append(sub[v].dropna().values)

    labels_units = {
        'wave_dir': ('Wave Direction', '[deg]'),
        'current_dir': ('Current Direction', '[deg]'),
        'water_level': ('Water Level', '[m]'),
        'Uc': ('Current Speed', '[m/s]'),
    }

    # Circular center for directional variables (puts this degree at mid-axis)
    dir_centers = {'wave_dir': 0, 'current_dir': 250}

    pal = sns.cubehelix_palette(n_top, rot=-.25, light=.7)
    rng = np.random.default_rng(42)

    for var in plot_vars:
        lname, unit = labels_units[var]
        datasets = bin_data[var]

        valid = [(i, d) for i, d in enumerate(datasets) if len(d) >= 2]
        if not valid:
            continue

        positions = [i for i, _ in valid]
        vals = [d for _, d in valid]

        # Circular shift for directional variables
        if var in dir_centers:
            c = dir_centers[var]
            vals = [(d - c + 180) % 360 - 180 + c for d in vals]

        fig, ax = plt.subplots(figsize=(10, 0.8 * n_top + 1.5))

        # half-violin
        parts = ax.violinplot(vals, positions=positions, vert=False,
                              showmedians=False, showextrema=False)
        for i, body in enumerate(parts['bodies']):
            m = np.mean(body.get_paths()[0].vertices[:, 1])
            body.get_paths()[0].vertices[:, 1] = np.clip(
                body.get_paths()[0].vertices[:, 1], m, np.inf)
            body.set_facecolor(pal[positions[i]])
            body.set_edgecolor('k')
            body.set_linewidth(0.5)
            body.set_alpha(0.8)

        # jittered scatter + stats (median, 95% CI)
        for i, (pos, d) in enumerate(zip(positions, vals)):
            jitter = rng.uniform(-0.3, 0.0, size=len(d))
            ax.scatter(d, pos + jitter, s=4, alpha=0.35, color=pal[pos],
                       edgecolors='none', zorder=3)

            med = np.median(d)
            ci_lo, ci_hi = np.percentile(d, [2.5, 97.5])
            ax.plot(med, pos, 'D', color='k', markersize=5, zorder=5)
            ax.hlines(pos, ci_lo, ci_hi, colors='k', linewidths=1.8, zorder=4)

        # x-axis sized to filtered data range with 5% margin
        all_vals = np.concatenate(vals)
        margin = (all_vals.max() - all_vals.min()) * 0.05
        ax.set_xlim(all_vals.min() - margin, all_vals.max() + margin)

        # Fix tick labels to show true degrees for directional variables
        if var in dir_centers:
            from matplotlib.ticker import FuncFormatter
            ax.xaxis.set_major_formatter(
                FuncFormatter(lambda x, _: f'{x % 360:.0f}'))

        ax.set_yticks(range(len(bin_labels)))
        ax.set_yticklabels(bin_labels, fontsize=8)
        ax.set_xlabel(f'{lname} {unit}', fontsize=11)
        ax.set_title(f'{lname} – Top {n_top} FLS bins (BME)', fontsize=12)
        ax.grid(axis='x', alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(graphs_dir, f'violin_{var}.png'),
                    dpi=150, bbox_inches='tight')
        plt.close()



if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    graphs_dir = os.path.join(script_dir, "graphs")
    os.makedirs(graphs_dir, exist_ok=True)

    metocean_file = os.path.join(script_dir, "metocean.xlsx")
    fls_file = os.path.join(script_dir, "fls_proxy_ff_v2.xlsx")

    hs_bins = 20
    tp_bins = 20
    n_sectors = 36

    df = load_data(metocean_file, fls_file)
    print(f"Loaded {len(df)} rows")
    print(df.head())

    plot_frequency(df, hs_bins, tp_bins, graphs_dir)
    plot_heatmaps(df, hs_bins, tp_bins, graphs_dir)
    plot_damage_roses(df, n_sectors, graphs_dir)

    plot_pareto_bins(df, hs_bins, tp_bins, n_top=30, graphs_dir=graphs_dir)
    print("Saved: ranked_fls_bins.png")

    plot_strain_n_curve(df, graphs_dir=graphs_dir)
    print("Saved: strain_n_curve.png")

    plot_violin_top_bins(df, metocean_file, hs_bins, tp_bins,
                         n_top=10, graphs_dir=graphs_dir)
    print("Saved: violin_*.png")
