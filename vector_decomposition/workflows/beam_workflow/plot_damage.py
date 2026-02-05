import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

MOMENT_SCALE = 1e-6   # divide moments by this (show in MN·m)
FLS_SCALE = 1e-12     # divide FLS by this

def load_data(metocean_path, fls_path):
    met_df = pd.read_excel(metocean_path, sheet_name='metocean')
    fls_df = pd.read_excel(fls_path)
    df = pd.merge(fls_df, met_df[['time', 'Hs', 'Tp', 'wave_dir']], on='time')
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
    # im = ax.imshow(H, origin='lower', aspect='auto',
    #                extent=[hs_edges[0], hs_edges[-1], tp_edges[0], tp_edges[-1]])
    im = ax.imshow(H, origin='lower', aspect='auto',
                   extent=[hs_edges[0], hs_edges[-1], tp_edges[0], tp_edges[-1]],
                   norm=PowerNorm(gamma=0.4), cmap='plasma')
    ax.set_xlabel('Hs [m]')
    ax.set_ylabel('Tp [s]')
    ax.set_title('Bin Frequency (count per Hs-Tp bin)')
    plt.colorbar(im, ax=ax, label='Count')
    plt.tight_layout()
    plt.savefig(os.path.join(graphs_dir, 'frequency.png'), dpi=150)

def plot_heatmaps(df, hs_bins=10, tp_bins=10, graphs_dir='graphs'):
    variables = [
        ('M_touchdown', f'Moment at TDP [×{MOMENT_SCALE:.0e}]', MOMENT_SCALE),
        ('M_bellmouth', f'Moment at BME [×{MOMENT_SCALE:.0e}]', MOMENT_SCALE),
        ('FLS_proxy_tdp', f'FLS Proxy at TDP [×{FLS_SCALE:.0e}]', FLS_SCALE),
        ('FLS_proxy_bme', f'FLS Proxy at BME [×{FLS_SCALE:.0e}]', FLS_SCALE)
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

        # im = ax.imshow(H, origin='lower', aspect='auto',
        #                extent=[hs_edges[0], hs_edges[-1], tp_edges[0], tp_edges[-1]])
        im = ax.imshow(H, origin='lower', aspect='auto',
                       extent=[hs_edges[0], hs_edges[-1], tp_edges[0], tp_edges[-1]],
                       norm=PowerNorm(gamma=0.4), cmap='plasma')
        ax.set_xlabel('Hs [m]')
        ax.set_ylabel('Tp [s]')
        ax.set_title(f'Cumulative {label}')
        plt.colorbar(im, ax=ax)

    plt.tight_layout()
    plt.savefig(os.path.join(graphs_dir, 'heatmaps.png'), dpi=150)
    # plt.show()

def plot_damage_roses(df, n_sectors=12, graphs_dir='graphs'):
    variables = [
        ('M_touchdown', f'Moment at TDP [×{MOMENT_SCALE:.0e}]', MOMENT_SCALE),
        ('M_bellmouth', f'Moment at BME [×{MOMENT_SCALE:.0e}]', MOMENT_SCALE),
        ('FLS_proxy_tdp', f'FLS Proxy at TDP [×{FLS_SCALE:.0e}]', FLS_SCALE),
        ('FLS_proxy_bme', f'FLS Proxy at BME [×{FLS_SCALE:.0e}]', FLS_SCALE)
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
    # plt.show()


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    graphs_dir = os.path.join(script_dir, "graphs")
    os.makedirs(graphs_dir, exist_ok=True)

    metocean_file = os.path.join(script_dir, "metocean.xlsx")
    fls_file = os.path.join(script_dir, "fls_proxy_global_local.xlsx")

    hs_bins = 50
    tp_bins = 50
    n_sectors = 36

    df = load_data(metocean_file, fls_file)
    print(f"Loaded {len(df)} rows")
    print(df.head())

    plot_frequency(df, hs_bins, tp_bins, graphs_dir)
    plot_heatmaps(df, hs_bins, tp_bins, graphs_dir)
    plot_damage_roses(df, n_sectors, graphs_dir)
