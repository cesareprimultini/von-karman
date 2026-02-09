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


def get_curvature_from_moment(M_demand):
    """Interpolates curvature from the Stick-Slip bending stiffness table."""
    M_bins = np.array([0, 1.676, 3.612, 8.918, 17.575]) * 1000  # Nm
    K_bins = np.array([0, 0.004, 0.016, 0.200, 0.500])           # 1/m
    return np.interp(abs(M_demand), M_bins, K_bins)


# ---------------------------------------------------------------------------
#  Force computation
# ---------------------------------------------------------------------------

def _segment_drag_force_vector(u_w, sign, wave_dir_vec, curr_vec,
                                seg_axis, seg_length, beam_axis,
                                rho, Cd, D):
    """Compute the cross-flow drag force on one segment as a 3-D VECTOR.

    Steps (same Morison physics as before, but direction is preserved):
        1. Total fluid velocity  = sign * u_w * wave_dir  +  current
        2. Decompose into components parallel and normal to the segment axis.
           Only the normal (cross-flow) component creates drag.
        3. Drag force  = 0.5 * rho * Cd * D * |V_n|^2 * L  in the direction of V_n
        4. Project the drag onto the plane perpendicular to the equivalent beam axis.

    Returns
    -------
    F_perp : ndarray, shape (3,)
        Drag force vector in the plane perpendicular to the beam axis.
    """
    # 1. total fluid velocity at segment midpoint
    V_total = sign * u_w * wave_dir_vec + curr_vec

    # 2. decompose relative to segment axis
    V_par  = np.dot(V_total, seg_axis) * seg_axis   # along pipe  (no drag)
    V_norm = V_total - V_par                         # cross-flow  (creates drag)
    V_n    = np.linalg.norm(V_norm)

    if V_n < 1e-12:
        return np.zeros(3)

    # 3. Morison cross-flow drag  (direction follows V_norm)
    F_drag = 0.5 * rho * Cd * D * V_n**2 * seg_length * (V_norm / V_n)

    # 4. project onto plane perpendicular to equivalent beam axis
    F_along = np.dot(F_drag, beam_axis) * beam_axis
    F_perp  = F_drag - F_along

    return F_perp


# ---------------------------------------------------------------------------
#  Beam solvers
# ---------------------------------------------------------------------------

def _solve_beam_unsigned(segments, u_w_all, idx, wave_dir_vec, curr_vec, sign,
                         beam_axis, load_positions, edge_positions, eval_points,
                         midspan, L_beam, rho, Cd, D):
    """Solve the beam with unsigned force magnitudes (legacy method).

    Used to obtain crest / trough absolute moments for reporting columns
    (M_bellmouth, M_touchdown, M_max, etc.).
    """
    loads = []
    for i, seg in enumerate(segments):
        u_w = u_w_all[seg['seg_id']][idx]
        F_perp = _segment_drag_force_vector(
            u_w, sign, wave_dir_vec, curr_vec,
            seg['axis'], seg['length'], beam_axis, rho, Cd, D)
        loads.append((load_positions[i], np.linalg.norm(F_perp)))

    M0, ML_val, R0 = fixed_fixed_beam(loads, L_beam)

    M_bme  = abs(M0)
    M_tdp  = abs(ML_val)
    M_mid  = abs(moment_at(midspan, loads, M0, R0))
    M_edges = [abs(moment_at(ep, loads, M0, R0)) for ep in edge_positions]
    M_max  = max(abs(moment_at(ep, loads, M0, R0)) for ep in eval_points)

    return {
        'M_bme': M_bme, 'M_tdp': M_tdp, 'M_mid': M_mid,
        'M_max': M_max, 'M_edges': M_edges
    }


def _solve_beam_force_range(segments, u_w_all, idx, wave_dir_vec, curr_vec,
                            beam_axis, load_positions, edge_positions,
                            eval_points, midspan, L_beam, rho, Cd, D):
    """Solve the beam once with FORCE RANGE loads (new method).

    For each segment the procedure is:
        F_crest  = drag force vector at crest  (sign = +1)
        F_trough = drag force vector at trough (sign = -1)
        F_range  = F_crest - F_trough          (vector subtraction)

    The vector subtraction correctly handles:
        - Force direction reversal when wave dominates current
        - Non-collinear wave and current directions
        - The case of zero current (pure wave oscillation)

    |F_range| is then used as the scalar beam load.  One beam solve gives
    the moment RANGE directly — no need to subtract two unsigned solves.
    """
    loads = []
    for i, seg in enumerate(segments):
        u_w = u_w_all[seg['seg_id']][idx]

        # drag force vector at wave crest (+u_w)
        F_crest = _segment_drag_force_vector(
            u_w, +1.0, wave_dir_vec, curr_vec,
            seg['axis'], seg['length'], beam_axis, rho, Cd, D)

        # drag force vector at wave trough (-u_w)
        F_trough = _segment_drag_force_vector(
            u_w, -1.0, wave_dir_vec, curr_vec,
            seg['axis'], seg['length'], beam_axis, rho, Cd, D)

        # force RANGE  =  crest vector  -  trough vector
        F_range_vec = F_crest - F_trough
        F_range_mag = np.linalg.norm(F_range_vec)

        loads.append((load_positions[i], F_range_mag))

    # one beam solve — moments ARE the moment range
    M0, ML_val, R0 = fixed_fixed_beam(loads, L_beam)

    M_range_bme = abs(M0)
    M_range_tdp = abs(ML_val)
    M_range_mid = abs(moment_at(midspan, loads, M0, R0))
    M_range_edges = [abs(moment_at(ep, loads, M0, R0)) for ep in edge_positions]
    M_range_max = max(abs(moment_at(ep, loads, M0, R0)) for ep in eval_points)

    return {
        'M_range_bme': M_range_bme, 'M_range_tdp': M_range_tdp,
        'M_range_mid': M_range_mid, 'M_range_max': M_range_max,
        'M_range_edges': M_range_edges
    }


# ---------------------------------------------------------------------------
#  Main computation
# ---------------------------------------------------------------------------

def compute_fls_proxy_ff(segments_path, velocities_path, output_path, cps_heading,
                          rho=1025, Cd=1.2, D=0.251, w_sub=500.0, EI_stick=419000.0):
    seg_df = pd.read_excel(segments_path)
    segments, bellmouth, touchdown = build_segments_3d(seg_df, cps_heading)

    n_seg = len(segments)
    L_arc = sum(seg['length'] for seg in segments)

    # Equivalent beam: straight-line chord from bellmouth to touchdown
    beam_vec = touchdown - bellmouth
    L_beam = np.linalg.norm(beam_vec)
    beam_axis = beam_vec / L_beam

    # Project segment midpoints onto beam axis (load application points)
    load_positions = []
    for seg in segments:
        a_i = np.dot(seg['midpoint'] - bellmouth, beam_axis)
        load_positions.append(max(0.0, min(a_i, L_beam)))

    # Project segment endpoints onto beam axis (edge positions)
    edge_positions = []
    for seg in segments:
        edge_positions.append(np.dot(seg['start'] - bellmouth, beam_axis))
    edge_positions.append(np.dot(segments[-1]['end'] - bellmouth, beam_axis))
    edge_positions = [max(0.0, min(ep, L_beam)) for ep in edge_positions]

    midspan = L_beam / 2

    # All candidate points for max|M|
    eval_points = sorted(set(edge_positions + load_positions + [midspan]))

    print(f"CPS heading: {cps_heading}")
    print(f"Bellmouth: {np.round(bellmouth, 2)}")
    print(f"Touchdown: {np.round(touchdown, 2)}")
    print(f"Segments: {n_seg}, Arc length: {L_arc:.2f} m, Beam chord: {L_beam:.2f} m\n")

    xl = pd.ExcelFile(velocities_path)
    sheet_names = [s for s in xl.sheet_names if s.startswith('seg_')]

    base_df = pd.read_excel(velocities_path, sheet_name=sheet_names[0])
    n_rows = len(base_df)

    u_w_all = {}
    for sheet in sheet_names:
        df = pd.read_excel(velocities_path, sheet_name=sheet)
        u_w_all[sheet] = df['u_w'].values

    # Check for current data
    has_current = 'Uc' in base_df.columns and not base_df['Uc'].isna().all()

    results = []
    for idx in range(n_rows):
        time_val = base_df.loc[idx, 'time']
        wave_dir = base_df.loc[idx, 'wave_dir']
        Tp = base_df.loc[idx, 'Tp']

        wave_dir_rad = np.radians(wave_dir)
        wave_dir_vec = np.array([np.sin(wave_dir_rad), np.cos(wave_dir_rad), 0])

        # Current vector
        curr_vec = np.zeros(3)
        if has_current:
            u_c = base_df.loc[idx, 'Uc']
            c_dir_rad = np.radians(base_df.loc[idx, 'current_dir'])
            curr_vec = u_c * np.array([np.sin(c_dir_rad), np.cos(c_dir_rad), 0])

        common_args = dict(
            segments=segments, u_w_all=u_w_all, idx=idx,
            wave_dir_vec=wave_dir_vec, curr_vec=curr_vec,
            beam_axis=beam_axis, load_positions=load_positions,
            edge_positions=edge_positions, eval_points=eval_points,
            midspan=midspan, L_beam=L_beam, rho=rho, Cd=Cd, D=D
        )

        # Legacy unsigned beam solves (for reporting crest/trough moments)
        crest  = _solve_beam_unsigned(sign=+1.0, **common_args)
        trough = _solve_beam_unsigned(sign=-1.0, **common_args)

        # Force-range beam solve (correct strain range)
        rng = _solve_beam_force_range(**common_args)

        # --- Strain calculations ---
        r_outer = D / 2.0

        # Crest / trough strains (from unsigned solves, for reporting)
        strain_crest_max  = get_curvature_from_moment(crest['M_max'])  * r_outer
        strain_trough_max = get_curvature_from_moment(trough['M_max']) * r_outer
        strain_crest_bme  = get_curvature_from_moment(crest['M_bme'])  * r_outer
        strain_trough_bme = get_curvature_from_moment(trough['M_bme']) * r_outer
        strain_crest_tdp  = get_curvature_from_moment(crest['M_tdp'])  * r_outer
        strain_trough_tdp = get_curvature_from_moment(trough['M_tdp']) * r_outer

        # Strain RANGE from force-range beam solve (the key improvement)
        delta_strain_max = get_curvature_from_moment(rng['M_range_max']) * r_outer
        delta_strain_bme = get_curvature_from_moment(rng['M_range_bme']) * r_outer
        delta_strain_tdp = get_curvature_from_moment(rng['M_range_tdp']) * r_outer

        # --- Damage calculation ---
        n_cycles = (3600.0 / Tp) if Tp > 0 else 0
        m = 3

        # Moment-based FLS (crest state)
        fls_tdp = n_cycles * crest['M_tdp']**m
        fls_bme = n_cycles * crest['M_bme']**m
        fls_max = n_cycles * crest['M_max']**m

        # Strain-based FLS (from force-range method)
        fls_strain_max = n_cycles * delta_strain_max**m
        fls_strain_bme = n_cycles * delta_strain_bme**m
        fls_strain_tdp = n_cycles * delta_strain_tdp**m

        row = {
            # Legacy columns (crest state values)
            'time': time_val,
            'M_touchdown': crest['M_tdp'],
            'M_bellmouth': crest['M_bme'],
            'FLS_proxy_tdp': fls_tdp,
            'FLS_proxy_bme': fls_bme,
            'M_midspan': crest['M_mid'],
            'M_max': crest['M_max'],
            'FLS_proxy_max': fls_max,
            # Crest/trough moment pairs
            'M_bellmouth_crest': crest['M_bme'],
            'M_bellmouth_trough': trough['M_bme'],
            'M_tdp_crest': crest['M_tdp'],
            'M_tdp_trough': trough['M_tdp'],
            'M_max_crest': crest['M_max'],
            'M_max_trough': trough['M_max'],
            # Strain outputs
            'Tp': Tp,
            'Strain_Crest': strain_crest_max,
            'Strain_Trough': strain_trough_max,
            'Strain_Range': delta_strain_max,
            'Strain_Crest_bme': strain_crest_bme,
            'Strain_Trough_bme': strain_trough_bme,
            'Strain_Range_bme': delta_strain_bme,
            'Strain_Crest_tdp': strain_crest_tdp,
            'Strain_Trough_tdp': strain_trough_tdp,
            'Strain_Range_tdp': delta_strain_tdp,
            # Strain-based FLS
            'FLS_Proxy': fls_strain_max,
            'FLS_Proxy_bme': fls_strain_bme,
            'FLS_Proxy_tdp': fls_strain_tdp,
        }
        for j, M_e in enumerate(crest['M_edges']):
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

    cps_heading = 65 # SUPER IMPORTANT!!!!!!

    segments_file = os.path.join(script_dir, "segments.xlsx")
    velocities_file = os.path.join(script_dir, "output.xlsx")
    output_file = os.path.join(script_dir, "fls_proxy_ff_v2.xlsx")

    compute_fls_proxy_ff(segments_file, velocities_file, output_file, cps_heading)
