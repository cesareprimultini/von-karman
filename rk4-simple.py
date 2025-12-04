"""
Simple viscous vortex method for Von Kármán vortex streets
NO boundary condition enforcement - just natural shedding and advection
Use this for calculating vortex wake velocity contributions at probe locations
Combine with separate potential flow and global flow calculations
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from tqdm import tqdm

# PARAMETERS
D = 4.88           # Cylinder diameter [m]
a = D / 2         # Radius
spacing = 40.0    # Center-to-center spacing [m]
U_inf = 1.0       # Freestream velocity [m/s]
Re = 100.0        # Reynolds number
St = 0.2          # Strouhal number
nu = U_inf * D / Re  # Kinematic viscosity

# Simulation parameters
dt = 0.03         # Time step
total_time = 300.0  # Total simulation time
shed_period = D / (St * U_inf)  # Shedding period T_shed

# Downstream vortex removal boundary
x_removal = 150.0  # Remove vortices beyond this x-coordinate

# Vortex parameters
sigma_0 = 0.1 * D
Gamma_mag = 2.0 * np.pi * St * U_inf * D  # |Gamma| with C_gamma = 2*pi
theta_sep = np.radians(80)   # Separation angle ~80 deg from rear

# Rotation and flow angles
rotation_angle = 30.0   # Cylinder cluster rotation [degrees]
flow_angle_metocean = 270.0  # Flow direction [degrees]
flow_angle = 270.0 - flow_angle_metocean

# Cylinders
s = spacing / 2.0
cylinders = [
    {'x': -s, 'y': +s},
    {'x': +s, 'y': +s},
    {'x': -s, 'y': -s},
    {'x': +s, 'y': -s},
]

# Measurement probe location
measure_point = (32.0, 5.0)

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
    x_local = cyl['x'] + a * np.cos(theta)
    y_local = cyl['y'] + a * np.sin(theta)
    x_shed, y_shed = apply_transforms(x_local, y_local, rotation_angle, flow_angle)

    return {'x': x_shed, 'y': y_shed, 'gamma': gamma, 'sigma': sigma_0, 'birth_t': t}

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
print("NO boundary enforcement - just natural shedding and advection")
num_steps = int(total_time / dt)

for step in tqdm(range(num_steps), desc="Simulating", mininterval=0.5, unit="step"):
    # 1. Shed new vortices
    for c_idx, cyl in enumerate(cylinders):
        if time >= next_shed_times[c_idx]:
            upper = (shed_counters[c_idx] % 2 == 0)
            new_vortex = shed_vortex(cyl, time, upper=upper)
            all_vortices.append(new_vortex)
            shed_counters[c_idx] += 1
            next_shed_times[c_idx] += shed_period
            shed_times.append(time)

    N = len(all_vortices)
    if N == 0:
        time += dt
        continue

    # 2. Update vortex core sizes (viscous diffusion)
    birth_times = np.array([v['birth_t'] for v in all_vortices])
    ages = time - birth_times
    sigmas = np.sqrt(sigma_0**2 + 4 * nu * ages)
    for i in range(N):
        all_vortices[i]['sigma'] = sigmas[i]

    # 3. Prepare arrays for RK4
    sources_x = np.array([v['x'] for v in all_vortices])
    sources_y = np.array([v['y'] for v in all_vortices])
    gammas = np.array([v['gamma'] for v in all_vortices])

    # 4. RK4 advection
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

    # 5. Remove vortices
    all_vortices = [v for v in all_vortices if v['x'] <= x_removal]
    all_vortices = remove_internal_vortices(cylinders, all_vortices, rotation_angle, flow_angle)

    # 6. Measure velocity at probe point
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
    circle = Circle((x_cyl, y_cyl), D/2,
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
ax.set_title(f'Simple Vortex Wake (D={D}m, Re={Re}, t={total_time}s)\n'
             f'NO boundary enforcement - for velocity amplification calculations')
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)
ax.legend(loc='upper right')

plt.tight_layout()
plt.savefig('4_cylinder_simple_wake.png', dpi=600, bbox_inches='tight')
print("\nPlot saved to '4_cylinder_simple_wake.png'")

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

print(f"Done! S/D = {spacing/D:.1f}")
print("\nThis vortex wake field can be combined with:")
print("  1. Potential flow (local cylinder effects)")
print("  2. Global flow (terrain/hill effects)")
print("  for total velocity amplification at probe locations")
