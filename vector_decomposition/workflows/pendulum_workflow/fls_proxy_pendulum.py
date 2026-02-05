import pandas as pd
import numpy as np
import os

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
        midpoint = (start + end) / 2.0
        axis = end - start
        axis = axis / np.linalg.norm(axis)
        segments.append({
            'seg_id': row['seg_id'],
            'length': row['length'],
            'start': start,
            'end': end,
            'midpoint': midpoint,
            'axis': axis
        })
    return segments

def get_curvature_from_moment(M_demand):
    """
    Interpolates curvature from the Stick-Slip table.
    Data Source: Image image_cdbf86.png
    """
    M_bins = np.array([0, 1.676, 3.612, 8.918, 17.575]) * 1000 
    K_bins = np.array([0, 0.004, 0.016, 0.200, 0.500])
    return np.interp(abs(M_demand), M_bins, K_bins)

def calculate_state_physics(segments, v_wave_vec, v_curr_vec, rho, Cd, D, W_total_sub, EI_stick):
    """
    Calculates the equilibrium angle, tension, and moment for a given velocity vector.
    """
    F_lat_vector = np.zeros(3)
    V_total = v_wave_vec + v_curr_vec # Vector addition (handles non-collinear)

    for seg in segments:
        # Cross-flow drag principle
        V_parallel = np.dot(V_total, seg['axis']) * seg['axis']
        V_normal = V_total - V_parallel
        V_n_mag = np.linalg.norm(V_normal)
        
        if V_n_mag > 0:
            F_seg_mag = 0.5 * rho * Cd * D * (V_n_mag**2) * seg['length']
            F_lat_vector += F_seg_mag * (V_normal / V_n_mag)

    F_lat_mag = np.linalg.norm(F_lat_vector)
    theta = np.arctan2(F_lat_mag, W_total_sub)
    T_top = W_total_sub / np.cos(theta)
    
    # Analytic Moment at bellmouth: M = theta * sqrt(T*EI)
    M_demand = theta * np.sqrt(T_top * EI_stick)
    
    return M_demand, theta, T_top

def compute_fls_proxy_range_logic(segments_path, velocities_path, output_path, cps_heading,
                                   rho=1025, Cd=1.2, D=0.251, w_sub=500.0):
    # 1. Setup
    seg_df = pd.read_excel(segments_path)
    segments = build_segments_3d(seg_df, cps_heading)
    total_length = sum(s['length'] for s in segments)
    W_total_sub = total_length * w_sub
    EI_stick = 419 * 1000  # Nm^2
    m_slope = 3.0 

    # 2. Load Metocean
    xl = pd.ExcelFile(velocities_path)
    sheet_names = [s for s in xl.sheet_names if s.startswith('seg_')]
    base_df = pd.read_excel(velocities_path, sheet_name=sheet_names[0])
    
    # Pre-load wave velocities
    u_w_all = {sheet: pd.read_excel(velocities_path, sheet_name=sheet)['u_w'].values for sheet in sheet_names}
    
    # Check for current data
    has_current = 'Uc' in base_df.columns and not base_df['Uc'].isna().all()
    
    n_rows = len(base_df)
    results = []

    print(f"Processing {n_rows} sea states using Range Logic (Once per cycle)...")

    for idx in range(n_rows):
        Tp = base_df.loc[idx, 'Tp']
        w_dir_rad = np.radians(base_df.loc[idx, 'wave_dir'])
        wave_unit_vec = np.array([np.sin(w_dir_rad), np.cos(w_dir_rad), 0])

        # Define Current Vector
        curr_vec = np.zeros(3)
        if has_current:
            u_c = base_df.loc[idx, 'Uc']
            c_dir_rad = np.radians(base_df.loc[idx, 'current_dir'])
            curr_vec = u_c * np.array([np.sin(c_dir_rad), np.cos(c_dir_rad), 0])

        # --- BOSS LOGIC: CALC TWO STATES PER WAVE ---
        # We assume u_w is the amplitude. Crest is +u_w, Trough is -u_w.
        u_w_amp = u_w_all[segments[0]['seg_id']][idx] # Using reference segment for simplicity
        wave_vec_crest = u_w_amp * wave_unit_vec
        wave_vec_trough = -u_w_amp * wave_unit_vec

        # Calculate Crest State (Max Flap)
        M_crest, theta_c, T_c = calculate_state_physics(segments, wave_vec_crest, curr_vec, rho, Cd, D, W_total_sub, EI_stick)
        kappa_crest = get_curvature_from_moment(M_crest)
        strain_crest = kappa_crest * (D / 2.0)

        # Calculate Trough State (Min Flap)
        M_trough, theta_t, T_t = calculate_state_physics(segments, wave_vec_trough, curr_vec, rho, Cd, D, W_total_sub, EI_stick)
        kappa_trough = get_curvature_from_moment(M_trough)
        strain_trough = kappa_trough * (D / 2.0)

        # --- DAMAGE CALCULATION ---
        # Delta Strain (Range)
        delta_strain = abs(strain_crest - strain_trough)
        
        # Count once per cycle (3600/Tp)
        n_cycles = (3600.0 / Tp) if Tp > 0 else 0
        damage_proxy = n_cycles * (delta_strain ** m_slope)

        results.append({
            'time': base_df.loc[idx, 'time'],
            'Tp': Tp,
            'Strain_Range': delta_strain,
            'Strain_Crest': strain_crest,
            'Strain_Trough': strain_trough,
            'Moment_Crest_kNm': M_crest / 1000,
            'Moment_Trough_kNm': M_trough / 1000,
            'FLS_Proxy': damage_proxy
        })

    # 3. Export
    output_df = pd.DataFrame(results)
    output_df.to_excel(output_path, index=False)
    print(f"Done. Range-based analysis saved to {output_path}")

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))

    cps_heading = 65

    segments_file = os.path.join(script_dir, "segments.xlsx")
    velocities_file = os.path.join(script_dir, "output.xlsx")
    output_file = os.path.join(script_dir, "fls_proxy_range_method.xlsx")

    compute_fls_proxy_range_logic(segments_file, velocities_file, output_file, cps_heading)