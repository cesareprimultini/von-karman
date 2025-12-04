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
dt = 0.03         # Time step (must be < sigma_0 / U_max ~0.5 for stability)
total_time = 300.0  # Total simulation time (adjust for longer wakes)
shed_period = D / (St * U_inf)  # Shedding period T_shed

# Downstream vortex removal boundary (to prevent unbounded accumulation)
x_removal = 150.0  # Remove vortices beyond this x-coordinate

# Vortex parameters 
sigma_0 = 0.1 * D
Gamma_mag = 2.0 * np.pi * St * U_inf * D  # |Gamma| with C_gamma = 2*pi
theta_sep = np.radians(80)   # Separation angle ~80 deg in radians from rear

# Rotation and flow angles
rotation_angle = 30.0   # Cylinder cluster rotation [degrees]
flow_angle_metocean = 270.0  # Flow direction [degrees, metocean: 0=from N, 90=from E, 180=from S, 270=from W]

# Convert metocean convention to mathematical convention (0=right, 90=up, 180=left, 270=down)
flow_angle = 270.0 - flow_angle_metocean

# Cylinders
s = spacing / 2.0
cylinders = [
    {'x': -s, 'y': +s},
    {'x': +s, 'y': +s},
    {'x': -s, 'y': -s},
    {'x': +s, 'y': -s},
]

# Apply rotation and flow angle transforms
def apply_transforms(x, y, rot_deg, flow_deg):
    total_angle = np.radians(rot_deg + flow_deg)
    c, s = np.cos(total_angle), np.sin(total_angle)
    return c*x - s*y, s*x + c*y

# Optimized function to compute velocities at arbitrary target points
def compute_velocities_at_points(targets_x, targets_y, sources_x, sources_y, gammas, sigmas):
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
    
    # Add contributions from all image vortices
    for cyl in cylinders:
        x_cyl, y_cyl = apply_transforms(cyl['x'], cyl['y'], rotation_angle, flow_angle)
        x_rel = sources_x - x_cyl
        y_rel = sources_y - y_cyl
        r_sq = x_rel**2 + y_rel**2
        mask_out = r_sq > a**2 + 1e-6
        x_img = np.full(N, 0.0)
        y_img = np.full(N, 0.0)
        x_img[mask_out] = x_cyl + (a**2 * x_rel[mask_out]) / r_sq[mask_out]
        y_img[mask_out] = y_cyl + (a**2 * y_rel[mask_out]) / r_sq[mask_out]
        gamma_img = -gammas
        sigma_img = sigmas
        
        for i in range(N):
            if not mask_out[i]:
                continue
            DX = targets_x - x_img[i]
            DY = targets_y - y_img[i]
            R = np.sqrt(DX**2 + DY**2)
            mask = R > 1e-6
            EXP = np.zeros_like(R)
            EXP[mask] = np.exp(-R[mask]**2 / sigma_img[i]**2)
            V_theta = np.zeros_like(R)
            V_theta[mask] = (gamma_img[i] / (2 * np.pi * R[mask])) * (1 - EXP[mask])
            U_ind = np.zeros_like(R)
            V_ind = np.zeros_like(R)
            U_ind[mask] = -V_theta[mask] * DY[mask] / R[mask]
            V_ind[mask] = V_theta[mask] * DX[mask] / R[mask]
            U += U_ind
            V += V_ind
    
    return U, V

# Shed vortex
def shed_vortex(cyl, t, upper=True):
    sign = 1 if upper else -1
    gamma = -Gamma_mag if upper else Gamma_mag  # upper clockwise (neg), lower counter (pos)

    theta = np.pi - sign * theta_sep
    x_local = cyl['x'] + a * np.cos(theta)
    y_local = cyl['y'] + a * np.sin(theta)
    x_shed, y_shed = apply_transforms(x_local, y_local, rotation_angle, flow_angle)

    return {'x': x_shed, 'y': y_shed, 'gamma': gamma, 'sigma': sigma_0, 'birth_t': t, 'parent_cyl': cyl}

# Main simulation (optimized)
all_vortices = []
time = 0.0
shed_counters = [0 for _ in cylinders]
next_shed_times = [0.0 for _ in cylinders]

print("Starting simulation...")
num_steps = int(total_time / dt)
for step in tqdm(range(num_steps), desc="Simulating", mininterval=0.5, unit="step"):
    # Shed new vortices
    for c_idx, cyl in enumerate(cylinders):
        if time >= next_shed_times[c_idx]:
            upper = (shed_counters[c_idx] % 2 == 0)
            new_vortex = shed_vortex(cyl, time, upper=upper)
            all_vortices.append(new_vortex)
            shed_counters[c_idx] += 1
            next_shed_times[c_idx] += shed_period
    
    N = len(all_vortices)
    if N == 0:
        time += dt
        continue
    
    # Vectorized sigma update
    birth_times = np.array([v['birth_t'] for v in all_vortices])
    ages = time - birth_times
    sigmas = np.sqrt(sigma_0**2 + 4 * nu * ages)
    for i in range(N):
        all_vortices[i]['sigma'] = sigmas[i]
    
    # Prepare arrays for sources (fixed during step)
    sources_x = np.array([v['x'] for v in all_vortices])
    sources_y = np.array([v['y'] for v in all_vortices])
    gammas = np.array([v['gamma'] for v in all_vortices])
    
    # RK4 advection (vectorized per stage)
    # Stage 1: at current positions
    U1, V1 = compute_velocities_at_points(sources_x, sources_y, sources_x, sources_y, gammas, sigmas)
    k1x = dt * U1
    k1y = dt * V1
    
    # Stage 2: at mid positions
    mid_x = sources_x + 0.5 * k1x
    mid_y = sources_y + 0.5 * k1y
    U2, V2 = compute_velocities_at_points(mid_x, mid_y, sources_x, sources_y, gammas, sigmas)
    k2x = dt * U2
    k2y = dt * V2
    
    # Stage 3
    mid_x = sources_x + 0.5 * k2x
    mid_y = sources_y + 0.5 * k2y
    U3, V3 = compute_velocities_at_points(mid_x, mid_y, sources_x, sources_y, gammas, sigmas)
    k3x = dt * U3
    k3y = dt * V3
    
    # Stage 4
    end_x = sources_x + k3x
    end_y = sources_y + k3y
    U4, V4 = compute_velocities_at_points(end_x, end_y, sources_x, sources_y, gammas, sigmas)
    k4x = dt * U4
    k4y = dt * V4
    
    # Update positions
    new_x = sources_x + (k1x + 2*k2x + 2*k3x + k4x) / 6
    new_y = sources_y + (k1y + 2*k2y + 2*k3y + k4y) / 6
    for i in range(N):
        all_vortices[i]['x'] = new_x[i]
        all_vortices[i]['y'] = new_y[i]

    # Remove vortices beyond downstream boundary
    all_vortices = [v for v in all_vortices if v['x'] <= x_removal]

    time += dt

print(f"Simulation complete. Total vortices: {len(all_vortices)}")

# Grid and velocity field computation (improved resolution)
grid_size = 600  
x_min = -80.0
x_max = 150.0  
y_min = -50.0
y_max = 50.0
x = np.linspace(x_min, x_max, grid_size)
y = np.linspace(y_min, y_max, grid_size)
X, Y = np.meshgrid(x, y)

def compute_velocity_field(X, Y, vortices):
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
        
        for cyl in cylinders:
            x_cyl, y_cyl = apply_transforms(cyl['x'], cyl['y'], rotation_angle, flow_angle)
            x_rel = v['x'] - x_cyl
            y_rel = v['y'] - y_cyl
            r_sq = x_rel**2 + y_rel**2
            if r_sq <= a**2 + 1e-6:
                continue
            x_img = x_cyl + (a**2 * x_rel) / r_sq
            y_img = y_cyl + (a**2 * y_rel) / r_sq
            gamma_img = -v['gamma']
            sigma_img = v['sigma']
            
            DX_img = X - x_img
            DY_img = Y - y_img
            R_img = np.sqrt(DX_img**2 + DY_img**2)
            mask_img = R_img > 1e-6
            EXP_img = np.zeros_like(R_img)
            EXP_img[mask_img] = np.exp(-R_img[mask_img]**2 / sigma_img**2)
            V_theta_img = np.zeros_like(R_img)
            V_theta_img[mask_img] = (gamma_img / (2 * np.pi * R_img[mask_img])) * (1 - EXP_img[mask_img])
            U_img = np.zeros_like(R_img)
            V_img = np.zeros_like(R_img)
            U_img[mask_img] = -V_theta_img[mask_img] * DY_img[mask_img] / R_img[mask_img]
            V_img[mask_img] = V_theta_img[mask_img] * DX_img[mask_img] / R_img[mask_img]
            U += U_img
            V += V_img
    
    return U, V

print("Computing velocity field...")
U, V = compute_velocity_field(X, Y, all_vortices)
vel_mag = np.sqrt(U**2 + V**2)

# PLOT (improved resolution and clarity)
fig, ax = plt.subplots(figsize=(16, 8), dpi=300)  # Larger figure with higher DPI

contour = ax.contourf(X, Y, vel_mag, levels=100, cmap='viridis')  # More levels for detail
plt.colorbar(contour, ax=ax, label='Velocity Magnitude [m/s]')

skip = 4  # Show more vectors with higher resolution grid
ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
          U[::skip, ::skip], V[::skip, ::skip],
          color='white', alpha=0.7, scale=100, width=0.0005)  # Smaller, finer arrows

for cyl in cylinders:
    x_cyl, y_cyl = apply_transforms(cyl['x'], cyl['y'], rotation_angle, flow_angle)
    circle = Circle((x_cyl, y_cyl), D/2,
                    color='gray', fill=True, zorder=10,
                    edgecolor='black', linewidth=2)
    ax.add_patch(circle)

for v in all_vortices:
    color = 'ro' if v['gamma'] > 0 else 'bo'
    ax.plot(v['x'], v['y'], color, markersize=2, alpha=0.4)  # Smaller markers for clarity

ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title(f'4-Cylinder Dynamic Wake (D={D}m, Spacing={spacing}m, Re={Re}, Rotation={rotation_angle}°, Flow={flow_angle_metocean}°, t={total_time}s)')
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('4_cylinder_dynamic_wake.png', dpi=600, bbox_inches='tight')  # High-res output
print("\nPlot saved to '4_cylinder_dynamic_wake.png' (high resolution)")
print(f"Done! S/D = {spacing/D:.1f}")