# sensitivity check on fls proxy outputs
#
# as per Terry email: we want to know which metocean input (hs, tp, wave dir, etc.)
# has the biggest effect on fatigue damage. instead of just correlation
# (which should miss nonlinear stuff), I chose to use a random forest classifier to predict
# fls from the inputs, then randomly shuffle each input one at a time.
# if shuffling a parameter lowers the prediction accuracy, that parameter
# matters a lot. thats bascally permutation importance.
# this is a bit of a black box method but I used this frequently in the past for sensitivity checks
# before training ML models.
#
# the R^2 on the plot title tells how well the model fits overall.
# if R^2 is low the ranking is unreliable. bars show mean importance
# across repeated shuffles, error bars show spread.
#
# two sets of features are used:
#   RAW_FEATURES  - columns straight from the metocean dataset
#   DERIVED_FEATURES - physically meaningful combinations computed from the raw data
#     angles are relative to the pipe axis (cps_heading=65 deg), not wave direction,
#     because cross-flow drag depends on velocity normal to the pipe.

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import StandardScaler

# --- feature definitions ------------------------------------------------

CPS_HEADING = 65  # pipe heading in degrees

RAW_FEATURES = ['Hs', 'Tp', 'wave_dir', 'current_dir', 'water_level', 'Uc']

DERIVED_FEATURES = [
    'Uw',                # max wave particle velocity across segments
    'wave_angle_to_pipe',  # angle between wave dir and pipe axis (0=inline, 90=crossflow)
    'curr_angle_to_pipe',  # angle between current dir and pipe axis
    'Uc_inline_pipe',      # current component along pipe axis
    'Uc_crossflow_pipe',   # current component normal to pipe axis
]

TARGETS = [
    ('FLS_Proxy', 'FLS Proxy (max strain)'),
    ('FLS_Proxy_bme', 'FLS Proxy (BME)'),
    ('FLS_Proxy_tdp', 'FLS Proxy (TDP)'),
]


def load_uw_from_output(output_path):
    """Read u_w from every segment sheet in output.xlsx and return max across segments per timestep."""
    xl = pd.ExcelFile(output_path)
    seg_sheets = [s for s in xl.sheet_names if s.startswith('seg_')]
    uw_max = None
    time_col = None
    for sheet in seg_sheets:
        sdf = pd.read_excel(output_path, sheet_name=sheet)
        if time_col is None:
            time_col = sdf['time']
        if uw_max is None:
            uw_max = sdf['u_w'].values.copy()
        else:
            uw_max = np.maximum(uw_max, sdf['u_w'].values)
    return pd.DataFrame({'time': time_col, 'Uw': uw_max})


def _angle_to_pipe(direction_deg):
    """Angle between a direction and the pipe axis, mapped to 0-90 (0=inline, 90=crossflow)."""
    diff = (direction_deg - CPS_HEADING + 180) % 360 - 180
    return np.abs(diff).clip(upper=90)


def add_derived_features(df):
    """Compute physically meaningful derived features relative to the pipe axis."""
    df = df.copy()

    # angles relative to pipe (0=inline with pipe, 90=perpendicular/crossflow)
    df['wave_angle_to_pipe'] = _angle_to_pipe(df['wave_dir'])
    df['curr_angle_to_pipe'] = _angle_to_pipe(df['current_dir'])

    # current components relative to pipe axis
    diff_rad = np.radians((df['current_dir'] - CPS_HEADING + 180) % 360 - 180)
    df['Uc_inline_pipe'] = df['Uc'] * np.cos(diff_rad)
    df['Uc_crossflow_pipe'] = df['Uc'] * np.abs(np.sin(diff_rad))

    return df


def _prepare_data(df, metocean_path, output_path):
    """Merge metocean + Uw, compute derived features, return clean df."""
    met_df = pd.read_excel(metocean_path, sheet_name='metocean')
    need = [c for c in RAW_FEATURES if c not in df.columns]
    if need:
        df = pd.merge(df, met_df[['time'] + need], on='time', how='left')

    if output_path is not None and 'Uw' not in df.columns:
        uw_df = load_uw_from_output(output_path)
        df = pd.merge(df, uw_df, on='time', how='left')

    df = add_derived_features(df)
    return df


def plot_sensitivity_importance(df, metocean_path, output_path=None, graphs_dir='graphs'):
    """Single combined importance plot with raw + derived features."""
    df = _prepare_data(df, metocean_path, output_path)
    targets = [(t, lab) for t, lab in TARGETS if t in df.columns]

    all_features = [f for f in RAW_FEATURES + DERIVED_FEATURES if f in df.columns]
    df = df.dropna(subset=all_features)
    X = StandardScaler().fit_transform(df[all_features].values)

    fig, axes = plt.subplots(1, len(targets), figsize=(6 * len(targets), 5), squeeze=False)

    for ax, (target, label) in zip(axes.flatten(), targets):
        y = df[target].values
        rf = RandomForestRegressor(n_estimators=200, max_depth=12, random_state=42, n_jobs=-1)
        rf.fit(X, y)
        r2 = rf.score(X, y)

        perm = permutation_importance(rf, X, y, n_repeats=10, random_state=42, n_jobs=-1)

        imp = perm.importances_mean
        order = np.argsort(imp)[::-1]
        ax.barh([all_features[i] for i in order],
                [imp[i] for i in order],
                xerr=[perm.importances_std[i] for i in order],
                color=sns.cubehelix_palette(len(all_features), rot=-.25, light=.7),
                edgecolor='k', linewidth=0.4)
        ax.set_xlabel('Permutation importance')
        ax.set_title(f'{label}  (R²={r2:.2f})')

    plt.tight_layout()
    plt.savefig(os.path.join(graphs_dir, 'sensitivity_importance.png'),
                dpi=150, bbox_inches='tight')
    plt.close()


def plot_cross_correlation(df, metocean_path, output_path=None, graphs_dir='graphs'):
    df = _prepare_data(df, metocean_path, output_path)
    all_features = [f for f in RAW_FEATURES + DERIVED_FEATURES if f in df.columns]
    df = df.dropna(subset=all_features)

    corr = df[all_features].corr()

    fig, ax = plt.subplots(figsize=(len(all_features) * 1.1 + 1,
                                    len(all_features) * 0.9 + 1))
    mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
    sns.heatmap(corr, annot=True, fmt='.2f', cmap='RdBu_r',
                center=0, vmin=-1, vmax=1, linewidths=0.5,
                mask=mask, square=True,
                xticklabels=all_features, yticklabels=all_features, ax=ax)
    ax.set_title('Cross-correlation of metocean inputs (raw + derived)')
    plt.tight_layout()
    plt.savefig(os.path.join(graphs_dir, 'sensitivity_cross_correlation.png'),
                dpi=150, bbox_inches='tight')
    plt.close()


def plot_input_distributions(df, metocean_path, output_path=None, graphs_dir='graphs'):
    df = _prepare_data(df, metocean_path, output_path)
    all_features = [f for f in RAW_FEATURES + DERIVED_FEATURES if f in df.columns]
    df = df.dropna(subset=all_features)

    fig, axes = plt.subplots(1, len(all_features),
                             figsize=(2.5 * len(all_features), 4), squeeze=False)

    for ax, feat in zip(axes.flatten(), all_features):
        sns.violinplot(y=df[feat], ax=ax, inner='quartile',
                       color=sns.cubehelix_palette(1, rot=-.25, light=.7)[0])
        ax.set_title(feat)
        ax.set_ylabel('')

    plt.suptitle('Distribution of metocean inputs (raw + derived)', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(graphs_dir, 'sensitivity_distributions.png'),
                dpi=150, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    from plot_damage import load_data

    script_dir = os.path.dirname(os.path.abspath(__file__))
    graphs_dir = os.path.join(script_dir, 'graphs')
    os.makedirs(graphs_dir, exist_ok=True)

    metocean_file = os.path.join(script_dir, 'metocean.xlsx')
    fls_file = os.path.join(script_dir, 'fls_proxy_ff_v2.xlsx')
    output_file = os.path.join(script_dir, 'output.xlsx')

    df = load_data(metocean_file, fls_file)

    plot_sensitivity_importance(df, metocean_file, output_file, graphs_dir=graphs_dir)
    print('Saved: sensitivity_importance.png')

    plot_cross_correlation(df, metocean_file, output_file, graphs_dir=graphs_dir)
    print('Saved: sensitivity_cross_correlation.png')

    plot_input_distributions(df, metocean_file, output_file, graphs_dir=graphs_dir)
    print('Saved: sensitivity_distributions.png')
