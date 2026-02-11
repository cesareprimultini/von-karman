import os
import pandas as pd
import numpy as np
from scipy.optimize import fsolve

g = 9.81


def solve_dispersion(T, d):
    """Solve the linear dispersion relation for wavenumber k and wavelength L."""
    omega = 2 * np.pi / T
    def dispersion(k):
        return omega**2 - g * k * np.tanh(k * d)
    # INITIAL GUESS: deep-water approx (swap with 0.01 or 0 to match Pat code)
    k0 = omega**2 / g
    # k0 = 0.01
    k = fsolve(dispersion, k0)[0]
    L = 2 * np.pi / k
    return k, L


def compute_wave_velocities(Hs, Tp, d, z_seabed):
    """Compute horizontal and vertical wave particle velocities (Airy wave theory)."""
    z = z_seabed - d
    k, L = solve_dispersion(Tp, d)
    u_w = (np.pi * Hs * np.cosh(k * (d + z))) / (Tp * np.sinh(k * d))
    w_w = (np.pi * Hs * np.sinh(k * (d + z))) / (Tp * np.sinh(k * d))
    return u_w, w_w, k, L


def compute_current_velocity(u_c_da, z_seabed, d):
    """Compute current velocity at segment midpoint using 1/7th power law."""
    u_c = (8/7) * (u_c_da) * (z_seabed / d)**(1/7)
    return u_c


def process_metocean(metocean_path, segments_path, output_path, d):
    """Process metocean data and compute velocities for each segment."""
    met_df = pd.read_excel(metocean_path, sheet_name='metocean')
    seg_df = pd.read_excel(segments_path)

    print(f"Config: d = {d}m")
    print(f"\nSegments (from {segments_path}):")
    print(seg_df[['seg_id', 'z_mid', 'length', 'theta_vert']].to_string(index=False))
    print()

    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        for _, seg in seg_df.iterrows():
            seg_id = seg['seg_id']
            z_seabed = seg['z_mid']

            results = []
            for idx, row in met_df.iterrows():
                d_inst = d + row['water_level']

                u_w, w_w, k, L = compute_wave_velocities(row['Hs'], row['Tp'], d_inst, z_seabed)
                u_c = compute_current_velocity(row['Uc_DA'], z_seabed, d_inst)

                results.append({
                    'k': k, 
                    'L': L,
                    'u_w': u_w, 
                    'w_w': w_w,
                    'u_c': u_c,
                })

                if (idx + 1) % 50000 == 0:
                    print(f"  {seg_id}: {idx + 1} / {len(met_df)} rows")

            res_df = pd.DataFrame(results)
            output = pd.DataFrame({
                'time': met_df['time'],
                'Tp': met_df['Tp'],
                'u_w': res_df['u_w'],
                'w_w': res_df['w_w'],
                'wave_dir': met_df['wave_dir'],
                'u_c': res_df['u_c'],
                'current_dir': met_df['current_dir'],
                'k': res_df['k'],
                'L': res_df['L'],
            })

            output.to_excel(writer, sheet_name=str(seg_id), index=False)
            print(f"Completed {seg_id} (z_mid = {z_seabed:.2f}m)")

    print(f"\nOutput saved to {output_path}")


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))

    d = 53.5
    metocean_file = os.path.join(script_dir, "metocean.xlsx")
    segments_file = os.path.join(script_dir, "segments.xlsx")
    output_file = os.path.join(script_dir, "output.xlsx")

    process_metocean(metocean_file, segments_file, output_file, d)
