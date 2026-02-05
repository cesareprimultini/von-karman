import pandas as pd
import numpy as np


def build_segments_3d(seg_df, cps_heading):
    heading_rad = np.radians(cps_heading)
    segments = []

    for _, row in seg_df.iterrows():
        start = np.array([
            row['x_start'] * np.sin(heading_rad),
            row['x_start'] * np.cos(heading_rad),
            row['z_start']
        ])
        end = np.array([
            row['x_end'] * np.sin(heading_rad),
            row['x_end'] * np.cos(heading_rad),
            row['z_end']
        ])
        midpoint = np.array([
            row['x_mid'] * np.sin(heading_rad),
            row['x_mid'] * np.cos(heading_rad),
            row['z_mid']
        ])
        axis = end - start
        axis = axis / np.linalg.norm(axis)

        segments.append({
            'seg_id': row['seg_id'],
            'length': row['length'],
            'theta_vert': row['theta_vert'],
            'start': start,
            'end': end,
            'midpoint': midpoint,
            'axis': axis
        })

    bellmouth = segments[0]['start']
    touchdown = segments[-1]['end']
    return segments, bellmouth, touchdown


def fixed_fixed_beam(loads, L):
    """Fixed-fixed end moments and left shear reaction for point loads [(a, P), ...]."""
    M0, ML = 0.0, 0.0
    for (a, P) in loads:
        b = L - a
        M0 += -P * a * b**2 / L**2
        ML += -P * a**2 * b / L**2
    tail = sum(P * (L - a) for (a, P) in loads)
    R0 = (ML - M0 + tail) / L
    return M0, ML, R0


def moment_at(x, loads, M0, R0):
    """M(x) = M0 + R0*x - sum P*(x-a) for loads with a <= x."""
    M = M0 + R0 * x
    for (a, P) in loads:
        if a <= x:
            M -= P * (x - a)
    return M


def compute_fls_proxy_ff(segments_path, velocities_path, output_path, cps_heading,
                          rho=1025, Cd=1.2, D=0.1):
    seg_df = pd.read_excel(segments_path)
    segments, bellmouth, touchdown = build_segments_3d(seg_df, cps_heading)

    n_seg = len(segments)
    L = sum(seg['length'] for seg in segments)

    # arc-length to each segment midpoint (where loads act)
    load_positions = []
    acc = 0.0
    for seg in segments:
        load_positions.append(acc + seg['length'] / 2)
        acc += seg['length']

    # arc-length to each segment edge (n_seg+1 values, 0 to L)
    edge_positions = [0.0]
    acc = 0.0
    for seg in segments:
        acc += seg['length']
        edge_positions.append(acc)

    midspan = L / 2

    # all candidate points for max|M|: edges + load positions + midspan
    eval_points = sorted(set(edge_positions + load_positions + [midspan]))

    print(f"CPS heading: {cps_heading}")
    print(f"Bellmouth: {np.round(bellmouth, 2)}")
    print(f"Touchdown: {np.round(touchdown, 2)}")
    print(f"Segments: {n_seg}, Arc length L: {L:.2f} m\n")

    xl = pd.ExcelFile(velocities_path)
    sheet_names = [s for s in xl.sheet_names if s.startswith('seg_')]

    base_df = pd.read_excel(velocities_path, sheet_name=sheet_names[0])
    n_rows = len(base_df)

    u_w_all = {}
    for sheet in sheet_names:
        df = pd.read_excel(velocities_path, sheet_name=sheet)
        u_w_all[sheet] = df['u_w'].values

    results = []
    for idx in range(n_rows):
        time_val = base_df.loc[idx, 'time']
        wave_dir = base_df.loc[idx, 'wave_dir']
        Tp = base_df.loc[idx, 'Tp']

        wave_dir_rad = np.radians(wave_dir)
        wave_dir_vec = np.array([np.sin(wave_dir_rad), np.cos(wave_dir_rad), 0])

        # build point loads for beam model: normal drag magnitude at each arc-length position
        loads = []
        for i, seg in enumerate(segments):
            u_w = u_w_all[seg['seg_id']][idx]
            V_wave = u_w * wave_dir_vec

            V_parallel = np.dot(V_wave, seg['axis']) * seg['axis']
            V_normal = V_wave - V_parallel
            V_n = np.linalg.norm(V_normal)

            if V_n > 0:
                F_mag = 0.5 * rho * Cd * D * V_n**2 * seg['length']
            else:
                F_mag = 0.0

            loads.append((load_positions[i], F_mag))

        # fixed-fixed beam solution
        M0, ML_val, R0 = fixed_fixed_beam(loads, L)

        # end moments
        M_bme = abs(M0)
        M_tdp = abs(ML_val)

        # midspan
        M_mid = abs(moment_at(midspan, loads, M0, R0))

        # segment edge moments
        M_edges = [abs(moment_at(ep, loads, M0, R0)) for ep in edge_positions]

        # max |M| over all candidate points (piecewise-linear extrema occur at load/edge positions)
        M_max = max(abs(moment_at(ep, loads, M0, R0)) for ep in eval_points)

        n_cycles = 2 * (3600 / Tp) if Tp > 0 else 0
        m = 3
        fls_tdp = n_cycles * M_tdp**m
        fls_bme = n_cycles * M_bme**m
        fls_max = n_cycles * M_max**m

        row = {
            'time': time_val,
            'M_touchdown': M_tdp,
            'M_bellmouth': M_bme,
            'FLS_proxy_tdp': fls_tdp,
            'FLS_proxy_bme': fls_bme,
            'M_midspan': M_mid,
            'M_max': M_max,
            'FLS_proxy_max': fls_max,
        }
        for j, M_e in enumerate(M_edges):
            row[f'M_edge_{j}'] = M_e

        results.append(row)

        if (idx + 1) % 50000 == 0:
            print(f"Processed {idx + 1} / {n_rows} rows")

    output_df = pd.DataFrame(results)
    output_df.to_excel(output_path, index=False)
    print(f"\nOutput saved to {output_path}")
    return output_df


if __name__ == "__main__":
    import os
    script_dir = os.path.dirname(os.path.abspath(__file__))

    cps_heading = 65

    segments_file = os.path.join(script_dir, "segments.xlsx")
    velocities_file = os.path.join(script_dir, "output.xlsx")
    output_file = os.path.join(script_dir, "fls_proxy_ff.xlsx")

    compute_fls_proxy_ff(segments_file, velocities_file, output_file, cps_heading)
