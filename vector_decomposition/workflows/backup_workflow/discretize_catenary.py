import pandas as pd
import numpy as np

def discretize_catenary(catenary_path, output_path, n_segments):
    df = pd.read_excel(catenary_path)
    df = df.drop_duplicates().sort_values("X' BM (m)").reset_index(drop=True)

    x = df["X' BM (m)"].values
    z = df["Z' NB (m)"].values

    x_min, x_max = x.min(), x.max()
    x_edges = np.linspace(x_min, x_max, n_segments + 1)

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

    print(f"\nDiscretized into {n_segments} segments:")
    print(seg_df.to_string(index=False))
    print(f"\nTouchdown at z = {segments[-1]['z_end']:.3f}m above seabed")
    print(f"Saved to {output_path}")
    return seg_df


if __name__ == "__main__":
    import os
    script_dir = os.path.dirname(os.path.abspath(__file__))

    n_segments = 5
    catenary_file = os.path.join(script_dir, "catenary_profile_data.xlsx")
    output_file = os.path.join(script_dir, "segments.xlsx")

    discretize_catenary(catenary_file, output_file, n_segments)
