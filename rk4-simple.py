"""
Viscous vortex method for Von Kármán vortex streets shed from one or more cylinders.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from tqdm import tqdm

# PARAMETERS
D = 4.88           # Cylinder diameter [m]
spacing = 40.0    # Center-to-center spacing [m]
U_inf = 1.0       # Freestream velocity [m/s]
Re = 1000000.0        # Reynolds number (100 = laminar, 6.7e6 = transcritical offshore)

# Turbulence modeling flags
enable_eddy_viscosity = (Re > 1000)  # Enable turbulent diffusion for Re > 1000
enable_stochastic_shedding = (Re > 1000)  # Enable random fluctuations for Re > 1000
enable_core_saturation = (Re > 1000)  # Limit max core size for Re > 1000
sigma_max_factor = 0.5  # Maximum core size = sigma_max_factor * D_ref (RECHECK THIS VALUE)

# Simulation parameters
dt = 0.01         # Time step
total_time = 300.0  # Total simulation time

# Downstream vortex removal boundary
x_removal = 250.0  # Remove vortices beyond this x-coordinate

# Rotation and flow angles
rotation_angle = 30.0   # Cylinder cluster rotation [degrees] (use 0.0 for one cylinder else fucks up)
flow_angle_metocean = 270.0  # Flow direction [degrees]
flow_angle = 270.0 - flow_angle_metocean

# Cylinders - easily configurable
s = spacing / 2.0
cylinders = [
    {'x': -s, 'y': +s, 'D': D},
    {'x': +s, 'y': +s, 'D': D},
    {'x': -s, 'y': -s, 'D': D},
    {'x': +s, 'y': -s, 'D': D},
]

# cylinders = [
#     {'x': 0.0, 'y': 0.0, 'D': D},
# ]

# Reference diameter and derived parameters
D_ref = max(cyl['D'] for cyl in cylinders)
nu = U_inf * D_ref / Re  # Kinematic viscosity

# Strouhal number (Re-dependent, Roshko 1961 + extensions)
if Re < 47:
    St = 0.0  # No vortex shedding
elif Re < 300:
    St = 0.212 * (1 - 21.2 / Re)  # Laminar regime
elif Re < 2e5:
    St = 0.20  # Subcritical regime
elif Re < 3.5e6:
    log_Re = np.log10(Re)
    St = 0.20 + 0.06 * (log_Re - 5.3)  # Critical/supercritical transition
else:
    St = 0.30  # Transcritical regime

shed_period = D_ref / (St * U_inf)  # Shedding period T_shed
sigma_0 = 0.1 * D_ref
Gamma_mag = 2.0 * np.pi * St * U_inf * D_ref  # |Gamma| with C_gamma = 2*pi
theta_sep = np.radians(80)   # Separation angle ~80 deg from rear

# Measurement probe location
measure_point = (32.0, 5.0)


def compute_eddy_viscosity_field(x, y, D, U_inf, Re):
    # Turbulent viscosity coefficient (calibrated from experiments)
    if Re < 1000:
        C_eddy = 0.01  # Minimal turbulent contribution
    elif Re < 1e5:
        # Gradual increase in subcritical regime
        C_eddy = 0.02 + 0.01 * np.log10(Re / 1000)
    else:
        # High-Re asymptote (supercritical/transcritical)
        C_eddy = 0.10

    # Normalized coordinates
    x_norm = x / D
    y_norm = y / D

    # Wake profile: Gaussian in transverse direction
    sigma_wake = 0.5 + 0.1 * x_norm  # Wake width grows linearly
    f_transverse = np.exp(-(y_norm**2) / (2 * sigma_wake**2))

    # Exponential decay in streamwise direction
    f_streamwise = np.exp(-x_norm / 20.0)  # Decays over ~20 diameters

    # Combined wake function
    f_wake = f_transverse * f_streamwise

    # Total effective viscosity (molecular + turbulent)
    nu_molecular = U_inf * D / Re
    nu_turbulent = C_eddy * U_inf * D * f_wake

    return nu_molecular + nu_turbulent


def compute_formation_time(D, U_inf, Re):
    # Base formation number
    F_star = 4.0

    # Reynolds correction for turbulent flows
    if Re > 1e5:
        F_star *= 0.8  # 20% faster formation

    T_form = F_star * (D / U_inf)
    return T_form


def circulation_growth_sigmoid(t_age, T_form):
    if t_age < 0:
        return 0.0

    # Normalized time
    tau = t_age / T_form

    # Sigmoid parameters
    k = 10.0  # Steepness

    # Raw sigmoid centered at tau = 0.5
    sigmoid_raw = 1.0 / (1.0 + np.exp(-k * (tau - 0.5)))

    # Normalize to [0, 1] range
    g_min = 1.0 / (1.0 + np.exp(k * 0.5))
    g_max = 1.0 / (1.0 + np.exp(-k * 0.5))
    g_normalized = (sigmoid_raw - g_min) / (g_max - g_min)

    return np.clip(g_normalized, 0.0, 1.0)


def shed_vortex_with_turbulence(cyl, t, upper, Re, U_inf, D_ref, Gamma_mag, sigma_0, rotation_angle, flow_angle, theta_sep, enable_stochastic=True):
    sign = 1 if upper else -1
    gamma_base = -Gamma_mag if upper else Gamma_mag

    # Add turbulent fluctuation to circulation (for high Re)
    if enable_stochastic and Re > 1000:
        # Fluctuation strength increases with log(Re)
        fluctuation_level = 0.10 * np.log10(Re / 1000)
        fluctuation_level = min(fluctuation_level, 0.20)  # Cap at 20%

        # Random Gaussian perturbation
        gamma = gamma_base * (1 + np.random.normal(0, fluctuation_level))
    else:
        gamma = gamma_base

    # Compute separation angle with optional randomness
    if enable_stochastic and Re > 1000:
        theta_sep_perturbed = theta_sep + np.random.normal(0, np.radians(5))
    else:
        theta_sep_perturbed = theta_sep

    # Shedding location on cylinder surface
    theta = np.pi - sign * theta_sep_perturbed
    a = cyl['D'] / 2
    x_local = cyl['x'] + a * np.cos(theta)
    y_local = cyl['y'] + a * np.sin(theta)

    # Apply rotation and flow angle transforms
    x_shed, y_shed = apply_transforms(x_local, y_local, rotation_angle, flow_angle)

    return {
        'x': x_shed,
        'y': y_shed,
        'gamma_target': gamma,      # Maximum circulation (final value)
        'gamma': 0.0,                # Current circulation (starts at zero)
        'sigma': sigma_0,
        'birth_t': t,
        'T_form': compute_formation_time(D_ref, U_inf, Re)  # Formation time scale
    }

def apply_transforms(x, y, rot_deg, flow_deg):
    """Apply rotation and flow angle transforms"""
    total_angle = np.radians(rot_deg + flow_deg)
    c, s = np.cos(total_angle), np.sin(total_angle)
    return c*x - s*y, s*x + c*y

def compute_velocities_at_points(targets_x, targets_y, sources_x, sources_y, gammas, sigmas):
    """Calculate velocities from vortices (NO image vortices, NO boundary enforcement)"""
    M = len(targets_x)
    if M == 0:
        return np.array([]), np.array([])

    U = np.full(M, U_inf * np.cos(np.radians(flow_angle)))
    V = np.full(M, U_inf * np.sin(np.radians(flow_angle)))
    N = len(sources_x)

    # Add contributions from all direct vortices
    for i in range(N):
        DX = targets_x - sources_x[i]
        DY = targets_y - sources_y[i]
        R = np.sqrt(DX**2 + DY**2)
        mask = R > 1e-6
        EXP = np.zeros_like(R)
        EXP[mask] = np.exp(-R[mask]**2 / sigmas[i]**2)
        V_theta = np.zeros_like(R)
        V_theta[mask] = (gammas[i] / (2 * np.pi * R[mask])) * (1 - EXP[mask])
        U_ind = np.zeros_like(R)
        V_ind = np.zeros_like(R)
        U_ind[mask] = -V_theta[mask] * DY[mask] / R[mask]
        V_ind[mask] = V_theta[mask] * DX[mask] / R[mask]
        U += U_ind
        V += V_ind

    return U, V

def shed_vortex(cyl, t, upper=True):
    """Shed a vortex from cylinder separation point"""
    sign = 1 if upper else -1
    gamma = -Gamma_mag if upper else Gamma_mag  # upper clockwise, lower counter

    theta = np.pi - sign * theta_sep
    a = cyl['D'] / 2
    x_local = cyl['x'] + a * np.cos(theta)
    y_local = cyl['y'] + a * np.sin(theta)
    x_shed, y_shed = apply_transforms(x_local, y_local, rotation_angle, flow_angle)

    return {
        'x': x_shed,
        'y': y_shed,
        'gamma_target': gamma,
        'gamma': 0.0,
        'sigma': sigma_0,
        'birth_t': t,
        'T_form': compute_formation_time(D_ref, U_inf, Re)
    }

def remove_internal_vortices(cylinders, vortices, rotation_angle, flow_angle):
    """Remove vortices inside cylinders"""
    filtered_vortices = []
    tolerance = 0.9  # Remove only if clearly inside

    for vortex in vortices:
        inside_any = False
        for cyl in cylinders:
            x_cyl, y_cyl = apply_transforms(cyl['x'], cyl['y'], rotation_angle, flow_angle)
            dx = vortex['x'] - x_cyl
            dy = vortex['y'] - y_cyl
            dist = np.sqrt(dx**2 + dy**2)
            a = cyl['D'] / 2
            if dist < a * tolerance:
                inside_any = True
                break
        if not inside_any:
            filtered_vortices.append(vortex)

    return filtered_vortices

# Main simulation
all_vortices = []
time = 0.0
shed_counters = [0 for _ in cylinders]
next_shed_times = [0.0 for _ in cylinders]
measure_times, measure_ux, measure_uy, measure_vmag = [], [], [], []
shed_times = []  # Track shedding times for visualization

print("Starting simple vortex method simulation...")
num_steps = int(total_time / dt)

for step in tqdm(range(num_steps), desc="Simulating", mininterval=0.5, unit="step"):
    # 1. Shed new vortices
    for c_idx, cyl in enumerate(cylinders):
        if time >= next_shed_times[c_idx]:
            upper = (shed_counters[c_idx] % 2 == 0)
            # Use turbulence-enabled shedding function
            new_vortex = shed_vortex_with_turbulence(
                cyl, time, upper, Re, U_inf, D_ref, Gamma_mag, sigma_0,
                rotation_angle, flow_angle, theta_sep, enable_stochastic=enable_stochastic_shedding
            )
            all_vortices.append(new_vortex)
            shed_counters[c_idx] += 1
            next_shed_times[c_idx] += shed_period
            shed_times.append(time)

    N = len(all_vortices)
    if N == 0:
        time += dt
        continue

    # 2. Update vortex core sizes (viscous diffusion with optional turbulence)
    birth_times = np.array([v['birth_t'] for v in all_vortices])
    ages = time - birth_times

    if enable_eddy_viscosity:
        # Compute spatially-varying effective viscosity for each vortex
        for i in range(N):
            x_vortex = all_vortices[i]['x']
            y_vortex = all_vortices[i]['y']
            nu_eff = compute_eddy_viscosity_field(x_vortex, y_vortex, D_ref, U_inf, Re)
            age = ages[i]
            sigma_new = np.sqrt(sigma_0**2 + 4 * nu_eff * age)

            # Apply core saturation if enabled
            if enable_core_saturation:
                sigma_max = sigma_max_factor * D_ref
                sigma_new = min(sigma_new, sigma_max)

            all_vortices[i]['sigma'] = sigma_new
    else:
        # Standard laminar diffusion (molecular viscosity only)
        sigmas_temp = np.sqrt(sigma_0**2 + 4 * nu * ages)
        for i in range(N):
            all_vortices[i]['sigma'] = sigmas_temp[i]

    # 3. Update vortex circulations (gradual formation)
    for i in range(N):
        t_age = time - all_vortices[i]['birth_t']
        growth_factor = circulation_growth_sigmoid(t_age, all_vortices[i]['T_form'])
        all_vortices[i]['gamma'] = all_vortices[i]['gamma_target'] * growth_factor

    # 4. Prepare arrays for RK4
    sources_x = np.array([v['x'] for v in all_vortices])
    sources_y = np.array([v['y'] for v in all_vortices])
    gammas = np.array([v['gamma'] for v in all_vortices])
    sigmas = np.array([v['sigma'] for v in all_vortices])  # Extract updated sigmas

    # 5. RK4 advection
    U1, V1 = compute_velocities_at_points(sources_x, sources_y, sources_x, sources_y, gammas, sigmas)
    k1x = dt * U1
    k1y = dt * V1

    mid_x = sources_x + 0.5 * k1x
    mid_y = sources_y + 0.5 * k1y
    U2, V2 = compute_velocities_at_points(mid_x, mid_y, sources_x, sources_y, gammas, sigmas)
    k2x = dt * U2
    k2y = dt * V2

    mid_x = sources_x + 0.5 * k2x
    mid_y = sources_y + 0.5 * k2y
    U3, V3 = compute_velocities_at_points(mid_x, mid_y, sources_x, sources_y, gammas, sigmas)
    k3x = dt * U3
    k3y = dt * V3

    end_x = sources_x + k3x
    end_y = sources_y + k3y
    U4, V4 = compute_velocities_at_points(end_x, end_y, sources_x, sources_y, gammas, sigmas)
    k4x = dt * U4
    k4y = dt * V4

    new_x = sources_x + (k1x + 2*k2x + 2*k3x + k4x) / 6
    new_y = sources_y + (k1y + 2*k2y + 2*k3y + k4y) / 6
    for i in range(N):
        all_vortices[i]['x'] = new_x[i]
        all_vortices[i]['y'] = new_y[i]

    # 6. Remove vortices
    all_vortices = [v for v in all_vortices if v['x'] <= x_removal]
    all_vortices = remove_internal_vortices(cylinders, all_vortices, rotation_angle, flow_angle)

    # 7. Measure velocity at probe point
    ux, uy = compute_velocities_at_points(
        np.array([measure_point[0]]),
        np.array([measure_point[1]]),
        sources_x, sources_y, gammas, sigmas
    )
    measure_times.append(time)
    measure_ux.append(ux[0] if len(ux) > 0 else U_inf * np.cos(np.radians(flow_angle)))
    measure_uy.append(uy[0] if len(uy) > 0 else U_inf * np.sin(np.radians(flow_angle)))
    measure_vmag.append(np.sqrt(measure_ux[-1]**2 + measure_uy[-1]**2))

    time += dt

print(f"Simulation complete. Total vortices in domain: {len(all_vortices)}")

# Compute velocity field
grid_size = 600
x_min = -80.0
x_max = 150.0
y_min = -50.0
y_max = 50.0
x = np.linspace(x_min, x_max, grid_size)
y = np.linspace(y_min, y_max, grid_size)
X, Y = np.meshgrid(x, y)

def compute_velocity_field(X, Y, vortices):
    """Compute velocity field from vortices only"""
    U = np.full_like(X, U_inf * np.cos(np.radians(flow_angle)))
    V = np.full_like(X, U_inf * np.sin(np.radians(flow_angle)))

    for v in vortices:
        DX = X - v['x']
        DY = Y - v['y']
        R = np.sqrt(DX**2 + DY**2)
        mask = R > 1e-6
        EXP = np.zeros_like(R)
        EXP[mask] = np.exp(-R[mask]**2 / v['sigma']**2)
        V_theta = np.zeros_like(R)
        V_theta[mask] = (v['gamma'] / (2 * np.pi * R[mask])) * (1 - EXP[mask])
        U_ind = np.zeros_like(R)
        V_ind = np.zeros_like(R)
        U_ind[mask] = -V_theta[mask] * DY[mask] / R[mask]
        V_ind[mask] = V_theta[mask] * DX[mask] / R[mask]
        U += U_ind
        V += V_ind

    return U, V

print("Computing velocity field...")
U, V = compute_velocity_field(X, Y, all_vortices)
vel_mag = np.sqrt(U**2 + V**2)

# Visualization
fig, ax = plt.subplots(figsize=(16, 8), dpi=300)

contour = ax.contourf(X, Y, vel_mag, levels=100, cmap='viridis')
plt.colorbar(contour, ax=ax, label='Velocity Magnitude [m/s]')

skip = 4
ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
          U[::skip, ::skip], V[::skip, ::skip],
          color='white', alpha=0.7, scale=100, width=0.0005)

for cyl in cylinders:
    x_cyl, y_cyl = apply_transforms(cyl['x'], cyl['y'], rotation_angle, flow_angle)
    circle = Circle((x_cyl, y_cyl), cyl['D']/2,
                    color='gray', fill=True, zorder=10,
                    edgecolor='black', linewidth=2)
    ax.add_patch(circle)

for v in all_vortices:
    color = 'ro' if v['gamma'] > 0 else 'bo'
    ax.plot(v['x'], v['y'], color, markersize=2, alpha=0.4)

ax.plot(measure_point[0], measure_point[1], 'r*', markersize=12,
        markeredgewidth=1.5, markeredgecolor='black', label='CPS', zorder=15)

ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title(f'Simple Vortex Wake (D={D_ref}m, Re={Re}, t={total_time}s)\n')
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)
ax.legend(loc='upper right')

plt.tight_layout()
plt.savefig('cylinder_wake.png', dpi=600, bbox_inches='tight')
print("\nPlot saved to 'cylinder_wake.png'")

# Velocity time history plot
print("Generating velocity measurement plot...")
fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(measure_times, measure_ux, label='$u_x$ (streamwise)',
        linewidth=1.5, alpha=0.8)
ax.plot(measure_times, measure_uy, label='$u_y$ (transverse)',
        linewidth=1.5, alpha=0.8)
ax.plot(measure_times, measure_vmag, label='$|v|$ (magnitude)',
        linewidth=2, color='black')

ax.set_xlabel('Time [s]', fontsize=12)
ax.set_ylabel('Velocity [m/s]', fontsize=12)
ax.set_title(f'Velocity Time History at ({measure_point[0]:.1f}, {measure_point[1]:.1f}) m',
             fontsize=13)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=11)
ax.set_xlim([0, total_time])
for t_shed in shed_times:
    ax.axvline(t_shed, color='red', alpha=0.15, linewidth=0.5, linestyle='-')

plt.tight_layout()
plt.savefig('velocity_measurement_simple.png', dpi=150, bbox_inches='tight')
print("Saved velocity_measurement_simple.png\n")

print(f"Done! S/D = {spacing/D_ref:.1f}")
print("\nThis vortex wake field can be combined with:")
print("  1. Potential flow (local cylinder effects)")
print("  2. Global flow (terrain/hill effects)")
print("  for total velocity amplification at probe locations")
