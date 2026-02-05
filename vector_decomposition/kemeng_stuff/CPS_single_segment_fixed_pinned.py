#!/usr/bin/env python3
"""
Fixed–Pinned (propped cantilever) equivalent beam under horizontal point loads
-----------------------------------------------------------------------------
Computes:
1) Left-end moment (fixed end) M0
2) Right-end moment ML = 0 (pinned)
3) Shear reactions at both ends (R0, RL)
4) Bending moment diagram M(x)
5) Candidate interior maxima/minima of M(x)

Model:
- Single equivalent Euler–Bernoulli beam
- x=0: fixed  -> w(0)=0, theta(0)=0 (reaction moment M0 allowed)
- x=L: pinned -> w(L)=0, M(L)=0 (rotation free)
- Horizontal point loads Pi at positions xi

Notes:
- This is still linear/small-deflection beam theory (a simple screening model).
"""

import numpy as np

# =========================
# INPUTS (edit these)
# =========================
L = 10.0  # [m] beam length

# Point loads: list of (x [m from left], P [N])
loads = [
    (5.0, 2000.0),
]

nx = 2001  # sampling resolution for moment diagram


# =========================
# CORE CALCULATIONS
# =========================
def M0_fixed_pinned_single_load(P, a, L):
    """
    Left fixed-end moment M0 for a fixed–pinned (propped cantilever) beam
    with a single point load P at x=a.

    Derivation via compatibility (M(L)=0 and w(L)=0) leads to:
        M0 = - P * a * (L - a) * (2L - a) / (2 L^2)

    Sanity checks:
    - a -> 0  => M0 -> 0
    - a -> L  => M0 -> 0
    """
    b = L - a
    return -P * a * b * (2*L - a) / (2 * L**2)


def M0_fixed_pinned(loads, L):
    """Superposition: total M0 is sum of individual contributions."""
    return sum(M0_fixed_pinned_single_load(P, a, L) for (a, P) in loads)


def end_reactions_from_ML(loads, L, M0, ML):
    """
    Given end moments M0 at x=0 and ML at x=L, compute end shears R0 and RL
    using the moment expression evaluated at x=L:

        M(L) = M0 + R0*L - sum(Pi*(L-xi)) = ML

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
    M(x) = M0 + R0*x - sum_{xi<=x} Pi*(x-xi)
    """
    x = np.asarray(x, dtype=float)
    M = M0 + R0 * x
    for (xi, P) in loads:
        dx = np.clip(x - xi, 0.0, None)
        M -= P * dx
    return M


def find_local_extrema(y):
    """Return indices of local extrema in sampled y."""
    dy = np.diff(y)
    eps = 1e-12 * max(1.0, np.max(np.abs(y)))
    s = np.sign(np.where(np.abs(dy) < eps, 0.0, dy))

    idx = []
    for i in range(1, len(s)):
        if s[i-1] > 0 and s[i] < 0:
            idx.append(i)
        elif s[i-1] < 0 and s[i] > 0:
            idx.append(i)
    return sorted(set(i for i in idx if 0 < i < len(y)-1))


# =========================
# RUN
# =========================
if __name__ == "__main__":
    # Validate
    for (x, P) in loads:
        if not (0.0 < x < L):
            raise ValueError(f"Load at x={x} must be within (0, L). Got L={L}.")

    # Fixed–Pinned boundary:
    ML = 0.0  # pinned end moment is zero by definition

    # 1) Compute left fixed-end moment M0 (superposition)
    M0 = M0_fixed_pinned(loads, L)

    # 2) Compute end reactions (shears)
    R0, RL = end_reactions_from_ML(loads, L, M0, ML)

    # 3) Moment diagram
    x = np.linspace(0.0, L, nx)
    M = moment_along_beam(x, loads, M0, R0)

    # 4) Candidates for interior peaks + endpoints
    cand = find_local_extrema(M)
    candidates = [(0.0, M[0], "end fixed (x=0)"), (L, M[-1], "end pinned (x=L)")]
    for i in cand:
        candidates.append((x[i], M[i], "interior"))
    candidates_sorted = sorted(candidates, key=lambda t: abs(t[1]), reverse=True)

    # =========================
    # PRINT
    # =========================
    print("\n=== Fixed–Pinned beam under horizontal point loads ===")
    print(f"L = {L:.6g} m")
    print("\nLoads (x [m], P [N]):")
    for (xi, P) in loads:
        print(f"  x = {xi:10.4f} m,  P = {P:12.4f} N")

    print("\n--- End moments ---")
    print(f"M_left  (fixed, x=0) : {M0: .6e} N·m")
    print(f"M_right (pinned,x=L) : {ML: .6e} N·m  (should be 0)")
    print(f"|M_left|             : {abs(M0): .6e} N·m")

    print("\n--- End reactions (shear) ---")
    print(f"R_left  at x=0 : {R0: .6e} N")
    print(f"R_right at x=L : {RL: .6e} N")
    print(f"Check sum V: R0+RL = {R0+RL:.6e} N,  sum(P) = {sum(P for _, P in loads):.6e} N")

    print("\n--- Boundary checks ---")
    print(f"M(0) = {M[0]: .6e} N·m  (should equal M_left)")
    print(f"M(L) = {M[-1]: .6e} N·m  (should be ~0 for pinned)")

    print("\n--- Largest moment candidates (by |M|) ---")
    topn = min(10, len(candidates_sorted))
    for k in range(topn):
        xi, Mi, tag = candidates_sorted[k]
        print(f"{k+1:2d}) x = {xi:10.4f} m,  M = {Mi: .6e} N·m,  |M| = {abs(Mi): .6e}  [{tag}]")

    # Save CSV
    out_csv = "moment_diagram_fixed_pinned.csv"
    np.savetxt(out_csv, np.column_stack([x, M]), delimiter=",", header="x_m,M_Nm", comments="")
    print(f"\nSaved moment diagram to: {out_csv}")

    # Plot (optional)
    try:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(x, M)
        plt.xlabel("x [m]")
        plt.ylabel("Bending moment M(x) [N·m]")
        plt.title("Moment diagram: fixed–pinned beam (horizontal loads)")
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    except Exception as e:
        print(f"\nPlot skipped: {e}")
