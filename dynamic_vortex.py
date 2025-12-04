"""
Simple von Kármán Vortex Street Simulator
Lamb-Oseen vortex model with time-dependent core growth
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# =============================================================================
# PARAMETERS
# =============================================================================
D = 1.0           # Cylinder diameter
U_inf = 1.5       # Freestream velocity
Re = 150.0        # Reynolds number
St = 0.2          # Strouhal number
nu = U_inf * D / Re  # Kinematic viscosity

sigma_0 = 0.1 * D    # Initial vortex core radius
Gamma_shed = 2 * np.pi * St * U_inf * D  # Circulation per vortex

# Shedding
T_shed = D / (St * U_inf)  # Period between vortices
theta_sep = 70 * np.pi/180  # Separation angle from rear (70 deg)

# Simulation
dt = 0.015         # Timestep
t_max = 30.0       # Total time (reduced for more frames in shorter time)
n_steps = int(t_max / dt)

# Output frames
n_frames = 40
frame_interval = max(1, n_steps // n_frames)

# Grid for velocity field
nx, ny = 500, 80
x_grid = np.linspace(-0.5*D, 12*D, nx)
y_grid = np.linspace(-2.5*D, 2.5*D, ny)
X, Y = np.meshgrid(x_grid, y_grid)

# =============================================================================
# LAMB-OSEEN VELOCITY FIELD
# =============================================================================
def lamb_oseen_velocity(x, y, x_vortex, y_vortex, gamma, sigma):
    dx = x - x_vortex
    dy = y - y_vortex
    r_sq = dx**2 + dy**2
    r = np.sqrt(r_sq + 1e-10)

    exp_factor = np.exp(-r_sq / (sigma**2 + 1e-10))
    vel_mag = (gamma / (2 * np.pi * r)) * (1 - exp_factor)

    u = -vel_mag * dy / r
    v = vel_mag * dx / r
    return u, v

# =============================================================================
# VORTEX ADVECTION
# =============================================================================
def get_velocity_at_point(x, y, vortices, t):
    u_total = U_inf
    v_total = 0.0

    for xv, yv, gamma, t_birth in vortices:
        dist = np.sqrt((x - xv)**2 + (y - yv)**2)
        if dist < 0.02 * D:  # Self-induction cutoff
            continue

        age = t - t_birth
        sigma = np.sqrt(sigma_0**2 + 4 * nu * age)

        u, v = lamb_oseen_velocity(x, y, xv, yv, gamma, sigma)
        u_total += u
        v_total += v

    return u_total, v_total

# =============================================================================
# MAIN SIMULATION
# =============================================================================
vortices = []  # List of (x, y, gamma, t_birth)
frames = []    # Saved velocity fields

shed_count = 0
t = 0.0

print(f"Simulation: {n_steps} steps, dt={dt:.4f}")
print(f"Shedding period: {T_shed:.3f}, Gamma: {Gamma_shed:.3f}\n")

for step in range(n_steps):
    # Shed new vortex when needed
    expected_count = int(t / T_shed)
    while shed_count < expected_count:
        is_upper = (shed_count % 2 == 0)

        # Position on cylinder surface at separation angle
        theta = np.pi - theta_sep if is_upper else np.pi + theta_sep
        x_shed = (D/2 + sigma_0) * np.cos(theta)  # Just outside cylinder surface
        y_shed = (D/2 + sigma_0) * np.sin(theta)

        # Circulation: upper=clockwise (negative), lower=CCW (positive)
        gamma = -Gamma_shed if is_upper else Gamma_shed

        vortices.append([x_shed, y_shed, gamma, t])
        shed_count += 1

        if shed_count <= 10:
            print(f"Shed vortex {shed_count}: {'UPPER' if is_upper else 'LOWER'} "
                  f"at ({x_shed:.2f}, {y_shed:.2f}), Gamma={gamma:+.2f}")

    # Advect existing vortices using simple Euler integration
    new_vortices = []
    for xv, yv, gamma, t_birth in vortices:
        u, v = get_velocity_at_point(xv, yv, vortices, t)

        xv_new = xv + u * dt
        yv_new = yv + v * dt

        # Keep if still in domain
        if xv_new < 15*D and abs(yv_new) < 4*D:
            new_vortices.append([xv_new, yv_new, gamma, t_birth])

    vortices = new_vortices

    # Save frame
    if step % frame_interval == 0 or step == n_steps - 1:
        U_field = np.full_like(X, U_inf)
        V_field = np.zeros_like(Y)

        for xv, yv, gamma, t_birth in vortices:
            age = t - t_birth
            sigma = np.sqrt(sigma_0**2 + 4 * nu * age)
            u, v = lamb_oseen_velocity(X, Y, xv, yv, gamma, sigma)
            U_field += u
            V_field += v

        frames.append({
            't': t,
            'U': U_field.copy(),
            'V': V_field.copy(),
            'vortices': [(xv, yv, gamma, t-t_birth) for xv, yv, gamma, t_birth in vortices]
        })

        if len(frames) % 10 == 0:
            print(f"Frame {len(frames)}: t={t:.2f}, N_vortices={len(vortices)}")

    t += dt

print(f"\nSimulation complete. {len(frames)} frames saved.\n")

# =============================================================================
# CALCULATE GLOBAL MIN/MAX FOR CONSISTENT COLORBAR
# =============================================================================
print("Computing global velocity bounds...")
vel_min = float('inf')
vel_max = float('-inf')

for frame in frames:
    vel_mag = np.sqrt(frame['U']**2 + frame['V']**2)
    vel_min = min(vel_min, vel_mag.min())
    vel_max = max(vel_max, vel_mag.max())

print(f"Global velocity range: [{vel_min:.3f}, {vel_max:.3f}]\n")

# =============================================================================
# SAVE FRAMES
# =============================================================================
os.makedirs('frames', exist_ok=True)

for i, frame in enumerate(frames):
    fig, ax = plt.subplots(figsize=(14, 6))

    # Velocity magnitude
    vel_mag = np.sqrt(frame['U']**2 + frame['V']**2)
    cf = ax.contourf(X, Y, vel_mag, levels=50, cmap='viridis', vmin=vel_min, vmax=vel_max)
    plt.colorbar(cf, ax=ax, label='Velocity |v| [m/s]')

    # Velocity vectors (higher resolution)
    skip = 4
    ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
              frame['U'][::skip, ::skip], frame['V'][::skip, ::skip],
              color='white', alpha=0.7, scale=50, width=0.0015, headwidth=3, headlength=4)

    # Vortex positions with growing cores
    for xv, yv, gamma, age in frame['vortices']:
        sigma = np.sqrt(sigma_0**2 + 4 * nu * age)
        color = 'red' if gamma < 0 else 'cyan'

        # Core radius circle
        circle = plt.Circle((xv, yv), sigma, color=color, fill=False,
                           linestyle=':', linewidth=1.2, alpha=0.7)
        ax.add_patch(circle)

        # Vortex center marker
        marker = '+' if gamma < 0 else 'x'
        ax.plot(xv, yv, marker, color=color, markersize=8, markeredgewidth=2)

    # Cylinder
    cylinder = plt.Circle((0, 0), D/2, facecolor='gray', edgecolor='black',
                         linewidth=2, zorder=10)
    ax.add_patch(cylinder)

    ax.set_xlim([x_grid[0], x_grid[-1]])
    ax.set_ylim([y_grid[0], y_grid[-1]])
    ax.set_xlabel('x/D', fontsize=11)
    ax.set_ylabel('y/D', fontsize=11)
    ax.set_title(f'Von Karman vortex street (Lamb-Oseen) - t={frame["t"]:.2f} s', fontsize=12)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2)

    filename = f'frames/frame_{i:03d}_t{frame["t"]:05.1f}.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()

print(f"Saved {len(frames)} frames to 'frames/' folder")
print(f"\nCreate video with:")
print(f"  ffmpeg -framerate 15 -i frames/frame_%03d_t*.png -c:v libx264 -pix_fmt yuv420p vortex_street.mp4")
