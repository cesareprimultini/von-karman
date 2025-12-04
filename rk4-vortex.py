import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from tqdm import tqdm

# PARAMETERS
D = 4.0           # Cylinder diameter [m]
a = D / 2         # Radius
spacing = 43.0    # Center-to-center spacing [m]
U_inf = 1.0       # Freestream velocity [m/s]
Re = 100.0        # Reynolds number
St = 0.2          # Strouhal number
nu = U_inf * D / Re  # Kinematic viscosity

# Simulation parameters
dt = 0.03         # Time step (must be < sigma_0 / U_max ~0.5 for stability)
total_time = 300.0  # Total simulation time (adjust for longer wakes)
shed_period = D / (St * U_inf)  # Shedding period T_shed

# Vortex parameters
sigma_0 = 0.1 * D
Gamma_mag = 2.0 * np.pi * St * U_inf * D  # |Gamma| with C_gamma = 2*pi
theta_sep = 1.4   # Separation angle ~80 deg in radians from rear

# Cylinders
s = spacing / 2
cylinders = [
    {'x': -s, 'y': +s},
    {'x': +s, 'y': +s},
    {'x': -s, 'y': -s},
    {'x': +s, 'y': -s},
]

# Function to compute image vortex
def get_image_vortex(vortex, cyl):
    x_rel = vortex['x'] - cyl['x']
    y_rel = vortex['y'] - cyl['y']
    r_sq = x_rel**2 + y_rel**2
    if r_sq <= a**2 + 1e-6:  # Slightly relaxed to avoid floating point issues on surface
        return None
    x_img_rel = (a**2 * x_rel) / r_sq
    y_img_rel = (a**2 * y_rel) / r_sq
    x_img = cyl['x'] + x_img_rel
    y_img = cyl['y'] + y_img_rel
    gamma_img = -vortex['gamma']
    sigma_img = vortex['sigma']
    return {'x': x_img, 'y': y_img, 'gamma': gamma_img, 'sigma': sigma_img}

# Scalar induced velocity
def induced_velocity(dx, dy, r, gamma, sigma):
    if r < 1e-6:
        return 0.0, 0.0
    exp_factor = np.exp(-r**2 / sigma**2)
    V_theta = (gamma / (2 * np.pi * r)) * (1 - exp_factor)
    u_ind = -V_theta * dy / r
    v_ind = V_theta * dx / r
    return u_ind, v_ind

# Total velocity at point (scalar, for advection)
def velocity_at_point(x, y, vortices, include_images=True, exclude_self_idx=-1):
    u = U_inf
    v = 0.0
    
    for idx, vortex in enumerate(vortices):
        if idx == exclude_self_idx:
            continue
        dx = x - vortex['x']
        dy = y - vortex['y']
        r = np.sqrt(dx**2 + dy**2)
        u_ind, v_ind = induced_velocity(dx, dy, r, vortex['gamma'], vortex['sigma'])
        u += u_ind
        v += v_ind
        
        if include_images:
            for cyl in cylinders:
                img = get_image_vortex(vortex, cyl)
                if img:
                    dx_img = x - img['x']
                    dy_img = y - img['y']
                    r_img = np.sqrt(dx_img**2 + dy_img**2)
                    u_img, v_img = induced_velocity(dx_img, dy_img, r_img, img['gamma'], img['sigma'])
                    u += u_img
                    v += v_img
    
    return u, v

# RK4 advection step
def rk4_advect(vortex_idx, vortices, dt):
    def vel_func(pos):
        return velocity_at_point(pos[0], pos[1], vortices, include_images=True, exclude_self_idx=vortex_idx)
    
    x0, y0 = vortices[vortex_idx]['x'], vortices[vortex_idx]['y']
    pos0 = np.array([x0, y0])
    
    k1 = np.array(vel_func(pos0)) * dt
    k2 = np.array(vel_func(pos0 + 0.5 * k1)) * dt
    k3 = np.array(vel_func(pos0 + 0.5 * k2)) * dt
    k4 = np.array(vel_func(pos0 + k3)) * dt
    
    pos_new = pos0 + (k1 + 2*k2 + 2*k3 + k4) / 6
    return pos_new[0], pos_new[1]

# Shed vortex
def shed_vortex(cyl, t, upper=True):
    sign = 1 if upper else -1
    gamma = -Gamma_mag if upper else Gamma_mag  # upper clockwise (neg), lower counter (pos)
    
    theta = np.pi - sign * theta_sep
    x_shed = cyl['x'] + a * np.cos(theta)
    y_shed = cyl['y'] + a * np.sin(theta)
    
    return {'x': x_shed, 'y': y_shed, 'gamma': gamma, 'sigma': sigma_0, 'birth_t': t, 'parent_cyl': cyl}

# Main simulation
all_vortices = []
time = 0.0
shed_counters = [0 for _ in cylinders]
next_shed_times = [0.0 for _ in cylinders]

print("Starting simulation...")
num_steps = int(total_time / dt)
with tqdm(total=num_steps, desc="Simulating", unit="step") as pbar:
    while time < total_time:
        # Shed
        for c_idx, cyl in enumerate(cylinders):
            if time >= next_shed_times[c_idx]:
                upper = (shed_counters[c_idx] % 2 == 0)
                new_vortex = shed_vortex(cyl, time, upper=upper)
                all_vortices.append(new_vortex)
                shed_counters[c_idx] += 1
                next_shed_times[c_idx] += shed_period

        # Update cores
        for v in all_vortices:
            age = time - v['birth_t']
            v['sigma'] = np.sqrt(sigma_0**2 + 4 * nu * age)

        # Advect
        new_positions = []
        for idx in range(len(all_vortices)):
            new_x, new_y = rk4_advect(idx, all_vortices, dt)
            new_positions.append((new_x, new_y))

        for idx, (new_x, new_y) in enumerate(new_positions):
            all_vortices[idx]['x'] = new_x
            all_vortices[idx]['y'] = new_y

        time += dt
        pbar.update(1)
        pbar.set_postfix({'vortices': len(all_vortices), 'time': f'{time:.1f}s'})

print(f"Simulation complete. Total vortices: {len(all_vortices)}")

# Grid
grid_size = 200
x_min = -50.0
x_max = 150.0  # Adjust based on total_time * U_inf
y_min = -50.0
y_max = 50.0
x = np.linspace(x_min, x_max, grid_size)
y = np.linspace(y_min, y_max, grid_size)
X, Y = np.meshgrid(x, y)

# Vectorized velocity field
def compute_velocity_field(X, Y, vortices):
    U = np.full_like(X, U_inf)
    V = np.full_like(X, 0.0)

    for v in tqdm(vortices, desc="Computing velocity field", unit="vortex"):
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
            x_rel = v['x'] - cyl['x']
            y_rel = v['y'] - cyl['y']
            r_sq = x_rel**2 + y_rel**2
            if r_sq <= a**2 + 1e-6:
                continue
            x_img = cyl['x'] + (a**2 * x_rel) / r_sq
            y_img = cyl['y'] + (a**2 * y_rel) / r_sq
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

# PLOT
fig, ax = plt.subplots(figsize=(12, 6))

contour = ax.contourf(X, Y, vel_mag, levels=50, cmap='viridis')
plt.colorbar(contour, ax=ax, label='Velocity Magnitude [m/s]')

skip = 4
ax.quiver(X[::skip, ::skip], Y[::skip, ::skip], 
          U[::skip, ::skip], V[::skip, ::skip],
          color='white', alpha=0.6, scale=20, width=0.002)

for cyl in cylinders:
    circle = Circle((cyl['x'], cyl['y']), D/2, 
                    color='gray', fill=True, zorder=10, 
                    edgecolor='black', linewidth=2)
    ax.add_patch(circle)

for v in all_vortices:
    color = 'ro' if v['gamma'] > 0 else 'bo'
    ax.plot(v['x'], v['y'], color, markersize=3, alpha=0.5)

ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title(f'4-Cylinder Dynamic Wake (D={D}m, Spacing={spacing}m, Re={Re}, t={total_time}s)')
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('4_cylinder_dynamic_wake.png')
print("\nPlot saved to '4_cylinder_dynamic_wake.png'")
print(f"Done! S/D = {spacing/D:.1f}")