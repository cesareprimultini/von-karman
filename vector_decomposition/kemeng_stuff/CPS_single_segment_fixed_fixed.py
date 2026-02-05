#!/usr/bin/env python3
"""
Fixed–fixed equivalent beam under horizontal point loads
-------------------------------------------------------
Computes:
1) End moments (fixed-end moments) at x=0 and x=L
2) Bending moment diagram M(x) along the beam
3) Locations and values of possible larger moments in the middle

Model:
- One equivalent Euler–Bernoulli beam, fixed at both ends
- Only horizontal point loads Pi applied at positions xi
- Uses closed-form fixed-end moment contributions + statics to build M(x)

Sign convention:
- Loads Pi > 0 act in +w direction (arbitrary).
- End moments returned are "reaction moments" (hogging negative by default in many texts).
- For "largest values", use absolute value |M|.

Edit the INPUTS section and run:
    python fixed_fixed_beam_pointloads.py
"""

import numpy as np

# =========================
# INPUTS (edit these)
# =========================
L = 10.0  # [m] equivalent beam length (bellmouth -> TDP)

# Point loads: list of (x [m from left end], P [N])
# Example: 3 loads along the beam
# loads = [
#     (2.0,  2000.0),   # (x, P)
#     (4.0, 3000.0),
#     (7.0, 1000.0),
# ]
loads = [
    (1,  2000.0),   # (x, P)
]

# Sampling resolution for the moment diagram
nx = 2001  # increase for smoother peak finding


# =========================
# CORE CALCULATIONS
# =========================
def fixed_end_moments_point_load(P, a, L):
    """
    Fixed-end moments for a fixed-fixed beam with a single point load P at x=a.
    b = L-a

    Returns (M_left, M_right) as reaction moments at x=0 and x=L.

    Common textbook form:
      M_left  = - P * a * b^2 / L^2
      M_right = - P * a^2 * b / L^2
    """
    b = L - a
    M_left = -P * a * b**2 / L**2
    M_right = -P * a**2 * b / L**2
    return M_left, M_right


def end_moments_all_loads(loads, L):
    """Superposition of fixed-end moments for multiple point loads."""
    ML, MR = 0.0, 0.0
    for (a, P) in loads:
        mL, mR = fixed_end_moments_point_load(P, a, L)
        ML += mL
        MR += mR
    return ML, MR


def end_reactions_from_compatibility(loads, L, M0, ML):
    """
    Compute end shears (reactions) for fixed-fixed beam given end moments M0 and ML
    and point loads, by enforcing M(L)=ML on the moment expression:

      M(x) = M0 + R0*x - sum_{xi <= x}(Pi*(x-xi))

    At x=L:
      ML = M0 + R0*L - sum(Pi*(L-xi))
      => R0 = (ML - M0 + sum(Pi*(L-xi))) / L
      RL = sum(Pi) - R0
    """
    Psum = sum(P for (x, P) in loads)
    tail = sum(P * (L - x) for (x, P) in loads)

    R0 = (ML - M0 + tail) / L
    RL = Psum - R0
    return R0, RL



def moment_along_beam(x, loads, M0, R0):
    """
    Bending moment M(x) for x in [0, L], using a left-section cut:

    M(x) = M0 + R0*x - sum_{xi <= x}(Pi*(x-xi))

    Where:
      - M0 is reaction moment at x=0
      - R0 is reaction shear at x=0
    """
    x = np.asarray(x, dtype=float)
    M = M0 + R0 * x
    for (xi, P) in loads:
        # contribution only for x >= xi
        dx = np.clip(x - xi, 0.0, None)
        M -= P * dx
    return M


def find_local_extrema(x, y):
    """
    Find local maxima/minima indices in y(x) by sign changes in first difference.
    Returns indices of candidates (excluding endpoints).
    """
    dy = np.diff(y)
    # Avoid tiny numerical noise
    eps = 1e-12 * max(1.0, np.max(np.abs(y)))
    s = np.sign(np.where(np.abs(dy) < eps, 0.0, dy))

    # Local extrema when slope changes sign: + to - (max), - to + (min)
    idx = []
    for i in range(1, len(s)):
        if s[i-1] > 0 and s[i] < 0:
            idx.append(i)      # i corresponds to y[i]
        elif s[i-1] < 0 and s[i] > 0:
            idx.append(i)
    # Remove endpoints just in case
    idx = [i for i in idx if 0 < i < len(y)-1]
    return sorted(set(idx))


# =========================
# RUN
# =========================
if __name__ == "__main__":
    # Validate load positions
    for (x, P) in loads:
        if not (0.0 < x < L):
            raise ValueError(f"Load at x={x} must be within (0, L). Got L={L}.")

    # 1) End moments from fixed-end formulas (superposition)
    M0, ML = end_moments_all_loads(loads, L)  # M0 at x=0, ML at x=L

    # 2) End shears from equilibrium
    R0, RL = end_reactions_from_compatibility(loads, L, M0, ML)


    # 3) Moment diagram
    x = np.linspace(0.0, L, nx)
    M = moment_along_beam(x, loads, M0, R0)

    # 4) Candidates for "larger values in the middle"
    #    For point loads, extrema can occur where shear crosses zero; our sampling finds local extrema.
    cand = find_local_extrema(x, M)

    # Always include endpoints as candidates too
    candidates = [(0.0, M[0], "end (x=0)"), (L, M[-1], "end (x=L)")]
    for i in cand:
        candidates.append((x[i], M[i], "interior"))

    # Sort by absolute moment magnitude descending
    candidates_sorted = sorted(candidates, key=lambda t: abs(t[1]), reverse=True)

    # =========================
    # PRINT RESULTS
    # =========================
    print("\n=== Fixed–fixed equivalent beam under horizontal point loads ===")
    print(f"L = {L:.6g} m")
    print("\nLoads (x [m], P [N]):")
    for (xi, P) in loads:
        print(f"  x = {xi:10.4f} m,  P = {P:12.4f} N")

    print("\n--- End moments (reaction moments) ---")
    print(f"M_left  at x=0 : {M0: .6e} N·m")
    print(f"M_right at x=L : {ML: .6e} N·m")
    print("\n(Use abs() if you only care about magnitude.)")
    print(f"|M_left|  : {abs(M0): .6e} N·m")
    print(f"|M_right| : {abs(ML): .6e} N·m")

    print("\n--- End reactions (shear, for consistency check) ---")
    print(f"R_left  at x=0 : {R0: .6e} N")
    print(f"R_right at x=L : {RL: .6e} N")
    print(f"Check sum V: R0+RL = {R0+RL:.6e} N,  sum(P) = {sum(P for _, P in loads):.6e} N")

    print("\n--- Largest moment candidates (by |M|) ---")
    topn = min(10, len(candidates_sorted))
    for k in range(topn):
        xi, Mi, tag = candidates_sorted[k]
        print(f"{k+1:2d}) x = {xi:10.4f} m,  M = {Mi: .6e} N·m,  |M| = {abs(Mi): .6e}  [{tag}]")

    # Also report the global max/min over the sampled grid
    i_max = int(np.argmax(M))
    i_min = int(np.argmin(M))
    print("\n--- Global extremes on sampled grid ---")
    print(f"Max  M: x = {x[i_max]:.4f} m,  M = {M[i_max]: .6e} N·m")
    print(f"Min  M: x = {x[i_min]:.4f} m,  M = {M[i_min]: .6e} N·m")
    print(f"Max |M|: x = {x[int(np.argmax(np.abs(M)))]:.4f} m,  |M| = {np.max(np.abs(M)):.6e} N·m")

    # Optional: write a CSV for plotting elsewhere
    out_csv = "moment_diagram_fixed_fixed.csv"
    data = np.column_stack([x, M])
    header = "x_m,M_Nm"
    np.savetxt(out_csv, data, delimiter=",", header=header, comments="")
    print(f"\nSaved moment diagram to: {out_csv}")

    # Optional quick plot (comment out if you don't want matplotlib)
    try:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(x, M)
        plt.xlabel("x [m]")
        plt.ylabel("Bending moment M(x) [N·m]")
        plt.title("Fixed–fixed beam moment diagram (horizontal loads)")
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    except Exception as e:
        print(f"\nPlot skipped (matplotlib not available or error): {e}")
