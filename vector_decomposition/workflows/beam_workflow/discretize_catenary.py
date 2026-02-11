import pandas as pd
import numpy as np


def _compute_curvature_weights(x, z):
    """Compute a weight at each point: 1 + |curvature|, normalized by arc length.

    Curved regions get higher weight, so segment edges cluster there.
    Flat regions still get weight >= 1 so they aren't starved of segments.
    """
    # Arc length between consecutive points
    ds = np.sqrt(np.diff(x)**2 + np.diff(z)**2)

    # Curvature via second derivative (finite differences)
    dz_dx = np.gradient(z, x)
    d2z_dx2 = np.gradient(dz_dx, x)
    curvature = np.abs(d2z_dx2) / (1 + dz_dx**2)**1.5

    # Weight: baseline 1 + curvature (so flat parts still get segments)
    w = 1.0 + curvature / (curvature.mean() + 1e-12)

    # Cumulative weighted arc length (trapezoidal integration)
    w_mid = (w[:-1] + w[1:]) / 2.0
    W = np.concatenate([[0], np.cumsum(w_mid * ds)])

    return W


def discretize_catenary(catenary_path, output_path, n_segments, tdp=None, adaptive=True):
    df = pd.read_excel(catenary_path)
    df = df.drop_duplicates().sort_values("X' BM (m)").reset_index(drop=True)

    # Filter to bellmouth-TDP region
    if tdp is not None:
        df = df[df["X' BM (m)"] <= tdp].reset_index(drop=True)
        print(f"Filtered catenary to TDP = {tdp}m ({len(df)} points retained)")

    # Shift profile up if any points are below seabed (e.g. scour hole) THIS NEEDS TO CHANGED SOON!!!
    z_min = df["Z' NB (m)"].min()
    if z_min < 0:
        df["Z' NB (m)"] = df["Z' NB (m)"] + abs(z_min)
        print(f"Shifted profile up by {abs(z_min):.3f}m (min Z was below seabed)")

    x = df["X' BM (m)"].values
    z = df["Z' NB (m)"].values

    x_min, x_max = x.min(), x.max()

    if adaptive and len(x) >= 4:
        # Place edges at equal intervals of curvature-weighted arc length
        W = _compute_curvature_weights(x, z)
        W_edges = np.linspace(0, W[-1], n_segments + 1)
        x_edges = np.interp(W_edges, W, x)
        # Force exact endpoints
        x_edges[0] = x_min
        x_edges[-1] = x_max
        mode = "adaptive"
    else:
        x_edges = np.linspace(x_min, x_max, n_segments + 1)
        mode = "uniform"

    segments = []
    for i in range(n_segments):
        x_start = x_edges[i]
        x_end = x_edges[i + 1]
        z_start = np.interp(x_start, x, z)
        z_end = np.interp(x_end, x, z)

        dx = x_end - x_start
        dz = z_start - z_end
        length = np.sqrt(dx**2 + dz**2)
        theta_vert = np.degrees(np.arctan2(dz, dx))

        segments.append({
            'seg_id': f'seg_{i+1}',
            'x_start': round(x_start, 4),
            'z_start': round(z_start, 4),
            'x_end': round(x_end, 4),
            'z_end': round(z_end, 4),
            'x_mid': round((x_start + x_end) / 2, 4),
            'z_mid': round((z_start + z_end) / 2, 4),
            'length': round(length, 4),
            'theta_vert': round(theta_vert, 2)
        })

    seg_df = pd.DataFrame(segments)
    seg_df.to_excel(output_path, index=False)

    print(f"\nDiscretized into {n_segments} segments ({mode}):")
    print(seg_df.to_string(index=False))
    print(f"\nTouchdown at z = {segments[-1]['z_end']:.3f}m above seabed")
    print(f"Saved to {output_path}")
    return seg_df


if __name__ == "__main__":
    import os
    script_dir = os.path.dirname(os.path.abspath(__file__))

    n_segments = 5
    tdp = 15  # horizontal distance from bellmouth to touchdown point (m)
    catenary_file = os.path.join(script_dir, "catenary_profile_data.xlsx")
    output_file = os.path.join(script_dir, "segments.xlsx")

    discretize_catenary(catenary_file, output_file, n_segments, tdp=tdp)
