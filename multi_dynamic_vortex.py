"""
Multi-Cylinder von Kármán Vortex Street Simulator
"""

import numpy as np
import matplotlib.pyplot as plt
import os

D = 4.88

cylinders = [
    {'x': 20, 'y': 20, 'D': D},
    {'x': 20, 'y': -20, 'D': D},
    {'x': -20, 'y': 20, 'D': D},
    {'x': -20, 'y': -20, 'D': D},
]

# cylinders = [
#     {'x': 0.0, 'y': 0.0, 'D': D},
# ]

rotation_angle = 30.0
flow_angle = 0.0

D_ref = max(cyl['D'] for cyl in cylinders)
U_inf = 1.5
Re = 150.0
St = 0.2
nu = U_inf * D_ref / Re

# Numerical regularization parameter scaled to characteristic length
# Reference: Cottet & Koumoutsakos, "Vortex Methods", Cambridge Univ Press, 2000
# Epsilon should be proportional to smallest physically relevant length scale
eps_sq = (1e-6 * D_ref)**2

sigma_0 = 0.1 * D_ref
Gamma_shed = 2 * np.pi * St * U_inf * D_ref
theta_sep = 70 * np.pi/180

measure_point = (32.0, 5.0)

dt = 0.005
t_max = 200.0
n_steps = int(t_max / dt)

n_frames = 40
frame_interval = max(1, n_steps // n_frames)

def lamb_oseen_velocity(x, y, x_vortex, y_vortex, gamma, sigma):
    dx = x - x_vortex
    dy = y - y_vortex
    r_sq = dx**2 + dy**2
    r = np.sqrt(r_sq + eps_sq)

    exp_factor = np.exp(-r_sq / (sigma**2 + eps_sq))
    vel_mag = (gamma / (2 * np.pi * r)) * (1 - exp_factor)

    u = -vel_mag * dy / r
    v = vel_mag * dx / r
    return u, v

def get_velocity_at_point(x, y, vortices, t):
    u_total = U_inf * np.cos(np.radians(flow_angle))
    v_total = U_inf * np.sin(np.radians(flow_angle))

    for xv, yv, gamma, t_birth, cyl_idx in vortices:
        age = t - t_birth
        sigma = np.sqrt(sigma_0**2 + 4 * nu * age)

        u, v = lamb_oseen_velocity(x, y, xv, yv, gamma, sigma)
        u_total += u
        v_total += v

    return u_total, v_total

def apply_transforms(x, y, rot_deg, flow_deg):
    total_angle = np.radians(rot_deg + flow_deg)
    c, s = np.cos(total_angle), np.sin(total_angle)
    return c*x - s*y, s*x + c*y

def compute_shed_position(cyl, is_upper, rot_deg, flow_deg):
    theta = np.pi - theta_sep if is_upper else np.pi + theta_sep
    x_local = cyl['x'] + (cyl['D']/2 + sigma_0) * np.cos(theta)
    y_local = cyl['y'] + (cyl['D']/2 + sigma_0) * np.sin(theta)
    return apply_transforms(x_local, y_local, rot_deg, flow_deg)

rotated_pos = [apply_transforms(c['x'], c['y'], rotation_angle, flow_angle) for c in cylinders]
xs, ys = zip(*rotated_pos)
Ds = [c['D'] for c in cylinders]

x_min = min(xs) - 2*max(Ds)
x_max = max(xs) + 15*max(Ds)
y_min = min(ys) - 3*max(Ds)
y_max = max(ys) + 3*max(Ds)

ny = 80
nx = int(ny * (x_max - x_min) / (y_max - y_min))
x_grid = np.linspace(x_min, x_max, nx)
y_grid = np.linspace(y_min, y_max, ny)
X, Y = np.meshgrid(x_grid, y_grid)

vortices = []
shed_counts = [0] * len(cylinders)
measure_times, measure_ux, measure_uy, measure_vmag = [], [], [], []
frames = []
shed_times = []  # check velocity bumps - delete

t = 0.0

print(f"Simulation: {n_steps} steps, dt={dt:.4f}")
print(f"Domain: x=[{x_min:.1f}, {x_max:.1f}], y=[{y_min:.1f}, {y_max:.1f}]")
print(f"Grid: {nx}x{ny}, {len(cylinders)} cylinders\n")

for step in range(n_steps):
    for i, cyl in enumerate(cylinders):
        T_shed = cyl['D'] / (St * U_inf)
        expected = int(t / T_shed)
        while shed_counts[i] < expected:
            is_upper = (shed_counts[i] % 2 == 0)
            x_shed, y_shed = compute_shed_position(cyl, is_upper, rotation_angle, flow_angle)
            gamma = -Gamma_shed if is_upper else Gamma_shed
            vortices.append([x_shed, y_shed, gamma, t, i])
            shed_times.append(t)  # check velocity bumps - delete
            shed_counts[i] += 1

            if shed_counts[i] <= 5:
                print(f"Cyl {i} shed {shed_counts[i]}: {'UPPER' if is_upper else 'LOWER'} "
                      f"at ({x_shed:.2f}, {y_shed:.2f}), Gamma={gamma:+.2f}")

    new_vortices = []
    for xv, yv, gamma, t_birth, cyl_idx in vortices:
        u, v = get_velocity_at_point(xv, yv, vortices, t)

        xv_new = xv + u * dt
        yv_new = yv + v * dt

        if xv_new < x_max and abs(yv_new) < abs(y_max):
            new_vortices.append([xv_new, yv_new, gamma, t_birth, cyl_idx])

    vortices = new_vortices

    ux, uy = get_velocity_at_point(measure_point[0], measure_point[1], vortices, t)
    measure_times.append(t)
    measure_ux.append(ux)
    measure_uy.append(uy)
    measure_vmag.append(np.sqrt(ux**2 + uy**2))

    if step % frame_interval == 0 or step == n_steps - 1:
        U_field = np.full_like(X, U_inf * np.cos(np.radians(flow_angle)))
        V_field = np.full_like(Y, U_inf * np.sin(np.radians(flow_angle)))

        for xv, yv, gamma, t_birth, cyl_idx in vortices:
            age = t - t_birth
            sigma = np.sqrt(sigma_0**2 + 4 * nu * age)
            u, v = lamb_oseen_velocity(X, Y, xv, yv, gamma, sigma)
            U_field += u
            V_field += v

        frames.append({
            't': t,
            'U': U_field.copy(),
            'V': V_field.copy(),
            'vortices': [(xv, yv, gamma, t-t_birth, cyl_idx) for xv, yv, gamma, t_birth, cyl_idx in vortices]
        })

        if len(frames) % 10 == 0:
            print(f"Frame {len(frames)}: t={t:.2f}, N_vortices={len(vortices)}")

    t += dt

print(f"\nSimulation complete. {len(frames)} frames saved.\n")

print("Computing global velocity bounds...")
vel_min = float('inf')
vel_max = float('-inf')

for frame in frames:
    vel_mag = np.sqrt(frame['U']**2 + frame['V']**2)
    vel_min = min(vel_min, vel_mag.min())
    vel_max = max(vel_max, vel_mag.max())

print(f"Global velocity range: [{vel_min:.3f}, {vel_max:.3f}]\n")

os.makedirs('frames', exist_ok=True)

skip = max(4, nx // 100)

for i, frame in enumerate(frames):
    fig, ax = plt.subplots(figsize=(14, 6))

    vel_mag = np.sqrt(frame['U']**2 + frame['V']**2)
    cf = ax.contourf(X, Y, vel_mag, levels=50, cmap='viridis', vmin=vel_min, vmax=vel_max)
    plt.colorbar(cf, ax=ax, label='Velocity |v| [m/s]')

    ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
              frame['U'][::skip, ::skip], frame['V'][::skip, ::skip],
              color='white', alpha=0.7, scale=50, width=0.0015, headwidth=3, headlength=4)

    for xv, yv, gamma, age, cyl_idx in frame['vortices']:
        sigma = np.sqrt(sigma_0**2 + 4 * nu * age)
        color = 'red' if gamma < 0 else 'cyan'

        circle = plt.Circle((xv, yv), sigma, color=color, fill=False,
                           linestyle=':', linewidth=1.2, alpha=0.7)
        ax.add_patch(circle)

        marker = '+' if gamma < 0 else 'x'
        ax.plot(xv, yv, marker, color=color, markersize=8, markeredgewidth=2)

    for cyl in cylinders:
        x_cyl, y_cyl = apply_transforms(cyl['x'], cyl['y'], rotation_angle, flow_angle)
        cylinder = plt.Circle((x_cyl, y_cyl), cyl['D']/2, facecolor='gray', edgecolor='black',
                             linewidth=2, zorder=10)
        ax.add_patch(cylinder)

    ax.plot(measure_point[0], measure_point[1], 'r*', markersize=12,
            markeredgewidth=1.5, markeredgecolor='black', label='CPS', zorder=15)

    ax.set_xlim([x_grid[0], x_grid[-1]])
    ax.set_ylim([y_grid[0], y_grid[-1]])
    ax.set_xlabel('x [m]', fontsize=11)
    ax.set_ylabel('y [m]', fontsize=11)
    ax.set_title(f'Multi-Cylinder von Karman - t={frame["t"]:.2f} s', fontsize=12)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2)
    ax.legend(loc='upper right')

    filename = f'frames/frame_{i:03d}_t{frame["t"]:05.1f}.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()

print(f"Saved {len(frames)} frames to 'frames/' folder\n")

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
ax.set_xlim([0, t_max])
for t_shed in shed_times: ax.axvline(t_shed, color='red', alpha=0.15, linewidth=0.5, linestyle='-')  # check velocity bumps - delete

plt.tight_layout()
plt.savefig('velocity_measurement.png', dpi=150, bbox_inches='tight')
print("Saved velocity_measurement.png\n")

# print(f"Create video with:")
# print(f"  ffmpeg -framerate 15 -i frames/frame_%03d_t*.png -c:v libx264 -pix_fmt yuv420p multi_vortex.mp4")
