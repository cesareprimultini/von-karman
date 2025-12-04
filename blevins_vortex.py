import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Parameters for the simulation
U = 1.0  # Free stream velocity (current from left to right)
a = 1.0  # Cylinder radius
sigma = 0.1  # Vortex core radius for desingularization
dt = 0.01  # Time step (reduced for accuracy)
t_max = 100.0  # Total simulation time (extended)
n_steps = int(t_max / dt) + 1

# Shedding parameters based on Strouhal number ~0.2 for cylinder
st = 0.2
d = 2 * a  # Diameter
f = st * U / d  # Shedding frequency
delta_t_shed = 1.0 / (2 * f)  # Time between consecutive sheds (alternate upper/lower)
shed_interval = int(delta_t_shed / dt)  # Steps between sheds
gamma_shed = 1.5  # Fixed shedding strength (arbitrary, adjust for realistic amplitude)
theta_sep = 1.4  # Separation angle in radians (~80 degrees from front stagnation)
delta_r = sigma  # Distance outside cylinder to place new vortex

# Precompute cylinder for plotting
theta_cyl = np.linspace(0, 2 * np.pi, 100)
x_cyl = a * np.cos(theta_cyl)
y_cyl = a * np.sin(theta_cyl)

# Grid for velocity field (extended and refined)
x_min, x_max = -2, 30
y_min, y_max = -5, 5
nx, ny = 100, 60  # Increased resolution
X, Y = np.meshgrid(np.linspace(x_min, x_max, nx), np.linspace(y_min, y_max, ny))
Z_grid = X + 1j * Y

# Save interval for animations (every 1 second)
save_interval = int(1.0 / dt)
velocity_fields = []  # List of (U_grid, V_grid, omega_grid)
field_times = []
positions_at_fields = []  # Current positions at save times

# Functions (vectorized)
def vortex_w(z, zeta, gamma, sigma):
    dz = z - zeta
    abs_dz = np.abs(dz)
    if np.isscalar(abs_dz):
        if abs_dz < 1e-8:
            return 0j
        else:
            return -1j * gamma * np.conj(dz) / (2 * np.pi * (sigma**2 + abs_dz**2))
    else:
        mask = abs_dz < 1e-8
        res = -1j * gamma * np.conj(dz) / (2 * np.pi * (sigma**2 + abs_dz**2))
        res[mask] = 0j
        return res

def uniform_w(z, U, a):
    abs_z = np.abs(z)
    return np.where(abs_z <= a, 0j, U * (1 - a**2 / z**2))

def compute_dz_dt(current_z, current_gammas, U, a, sigma):
    N = len(current_z)
    dz_dt = np.zeros(N, dtype=complex)
    for k in range(N):
        w_k = uniform_w(current_z[k], U, a)
        for j in range(N):
            if j == k:
                continue
            w_k += vortex_w(current_z[k], current_z[j], current_gammas[j], sigma)
        for j in range(N):
            z_img = a**2 / np.conj(current_z[j])
            w_k += vortex_w(current_z[k], z_img, -current_gammas[j], sigma)
        dz_dt[k] = np.conj(w_k)
    return dz_dt

# Initialization
current_z = []  # List of complex positions
current_gammas = []  # List of strengths
trajectories = []  # List of position histories for each vortex
birth_steps = []  # Step when each vortex was shed

# Simulation loop
for step in range(n_steps):
    t = step * dt

    # Shed new vortex periodically, alternating upper and lower
    if shed_interval > 0 and step % shed_interval == 0 and step > 0:
        is_upper = ((step // shed_interval) % 2 == 0)
        theta = np.pi - theta_sep if is_upper else np.pi + theta_sep
        sign = -1.0 if is_upper else 1.0  # Negative for upper (clockwise), positive for lower (CCW)
        z_new = (a + delta_r) * (np.cos(theta) + 1j * np.sin(theta))
        gamma_new = sign * gamma_shed
        current_z.append(z_new)
        current_gammas.append(gamma_new)
        trajectories.append([z_new])  # Start trajectory for new vortex
        birth_steps.append(step)

    # Update existing vortices using RK4
    N = len(current_z)
    if N > 0:
        current_z = np.array(current_z)
        current_gammas = np.array(current_gammas)
        
        k1 = dt * compute_dz_dt(current_z, current_gammas, U, a, sigma)
        k2 = dt * compute_dz_dt(current_z + k1 / 2, current_gammas, U, a, sigma)
        k3 = dt * compute_dz_dt(current_z + k2 / 2, current_gammas, U, a, sigma)
        k4 = dt * compute_dz_dt(current_z + k3, current_gammas, U, a, sigma)
        
        current_z += (k1 + 2 * k2 + 2 * k3 + k4) / 6
        
        # Append new positions to trajectories (for all vortices, including older ones)
        for k in range(N):
            trajectories[k].append(current_z[k])
        
        # Remove far-away vortices to save computation
        mask = np.real(current_z) < x_max + 5
        current_z = current_z[mask].tolist()
        current_gammas = current_gammas[mask].tolist()
        trajectories = [trajectories[i] for i in range(len(trajectories)) if mask[i]]
        birth_steps = [birth_steps[i] for i in range(len(birth_steps)) if mask[i]]

    # Save fields every 1 second
    if save_interval > 0 and step % save_interval == 0 and step > 0:
        W_grid = uniform_w(Z_grid, U, a)
        for j in range(N):
            W_grid += vortex_w(Z_grid, current_z[j], current_gammas[j], sigma)
            z_img = a**2 / np.conj(current_z[j])
            W_grid += vortex_w(Z_grid, z_img, -current_gammas[j], sigma)
        mask_inside = np.abs(Z_grid) <= a
        W_grid[mask_inside] = 0j
        U_grid = np.real(W_grid)
        V_grid = -np.imag(W_grid)
        
        # Compute vorticity
        dx = X[0, 1] - X[0, 0]
        dy = Y[1, 0] - Y[0, 0]
        dudy = np.gradient(U_grid, axis=0) / dy
        dvdx = np.gradient(V_grid, axis=1) / dx
        omega_grid = dvdx - dudy
        
        velocity_fields.append((U_grid.copy(), V_grid.copy(), omega_grid.copy()))
        field_times.append(t)
        positions_at_fields.append(current_z[:])

# Static plot: Vortex trajectories with wake boundary envelope
plt.figure(figsize=(12, 6))
plt.plot(x_cyl, y_cyl, 'k-', linewidth=2, label='Cylinder')
for traj in trajectories:
    traj = np.array(traj)
    plt.plot(np.real(traj), np.imag(traj), linewidth=1)

# Compute and plot wake boundary envelope
all_x = []
all_y = []
for traj in trajectories:
    traj = np.array(traj)
    all_x.extend(np.real(traj))
    all_y.extend(np.imag(traj))
all_x = np.array(all_x)
all_y = np.array(all_y)
if len(all_x) > 0:
    sort_idx = np.argsort(all_x)
    all_x = all_x[sort_idx]
    all_y = all_y[sort_idx]
    bins = np.linspace(np.min(all_x), np.max(all_x), 50)
    upper = []
    lower = []
    bin_centers = []
    for i in range(len(bins) - 1):
        mask = (all_x >= bins[i]) & (all_x < bins[i + 1])
        if np.sum(mask) > 0:
            ys = all_y[mask]
            upper.append(np.max(ys))
            lower.append(np.min(ys))
            bin_centers.append((bins[i] + bins[i + 1]) / 2)
    plt.plot(bin_centers, upper, 'r--', linewidth=2, label='Wake boundary upper')
    plt.plot(bin_centers, lower, 'r--', linewidth=2, label='Wake boundary lower')

plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Vortex Trajectories and Wake Boundary Envelope Over Time')
plt.legend()
plt.grid(True)
plt.show()

# Prepare for animations
num_frames = len(field_times)
save_steps = [int(t / dt) for t in field_times]

# Find global limits for consistent coloring
max_mag = max(np.max(np.sqrt(u**2 + v**2)) for u, v, _ in velocity_fields)
max_omega = max(np.max(np.abs(o)) for _, _, o in velocity_fields)

# Animation 1: Trajectory buildup over time
fig_traj, ax_traj = plt.subplots(figsize=(12, 6))
ax_traj.plot(x_cyl, y_cyl, 'k-', linewidth=2)
ax_traj.set_xlim(x_min, x_max)
ax_traj.set_ylim(y_min, y_max)
ax_traj.set_xlabel('x')
ax_traj.set_ylabel('y')
ax_traj.set_title('Vortex Trajectories Over Time')
ax_traj.grid(True)

lines = [ax_traj.plot([], [], linewidth=1)[0] for _ in trajectories]
scat_traj = ax_traj.scatter([], [], c='r', s=20, label='Current Vortices')

def init_traj():
    for line in lines:
        line.set_data([], [])
    scat_traj.set_offsets(np.empty((0, 2)))
    return lines + [scat_traj]

def update_traj(frame):
    current_step = save_steps[frame]
    current_pos = []
    for l, traj in enumerate(trajectories):
        birth = birth_steps[l]
        if current_step < birth:
            lines[l].set_data([], [])
            continue
        num_points = (current_step - birth) + 1  # +1 for initial
        num_points = min(num_points, len(traj))  # Safety
        x = np.real(traj[:num_points])
        y = np.imag(traj[:num_points])
        lines[l].set_data(x, y)
        current_pos.append(traj[num_points - 1])
    if current_pos:
        scat_traj.set_offsets(np.c_[np.real(current_pos), np.imag(current_pos)])
    ax_traj.set_title(f'Vortex Trajectories at t ≈ {field_times[frame]:.1f} s')
    return lines + [scat_traj]

ani_traj = FuncAnimation(fig_traj, update_traj, frames=num_frames, init_func=init_traj, blit=True)
ani_traj.save('trajectory_animation.mp4', writer='ffmpeg', fps=10)
plt.show()

# Animation 2: Velocity field (quiver) over time
fig_vel, ax_vel = plt.subplots(figsize=(12, 6))
ax_vel.plot(x_cyl, y_cyl, 'k-', linewidth=2)
ax_vel.set_xlim(x_min, x_max)
ax_vel.set_ylim(y_min, y_max)
ax_vel.set_xlabel('x')
ax_vel.set_ylabel('y')
ax_vel.set_title('Velocity Field Over Time')
ax_vel.grid(True)

# Initial quiver
u0, v0, _ = velocity_fields[0]
quiv = ax_vel.quiver(X, Y, u0, v0, scale=50, color='b')
scat_vel = ax_vel.scatter([], [], c='r', s=20, label='Vortices')

def init_vel():
    quiv.set_UVC(np.zeros_like(u0), np.zeros_like(v0))
    scat_vel.set_offsets(np.empty((0, 2)))
    return [quiv, scat_vel]

def update_vel(frame):
    u, v, _ = velocity_fields[frame]
    quiv.set_UVC(u, v)
    curr_pos = positions_at_fields[frame]
    scat_vel.set_offsets(np.c_[np.real(curr_pos), np.imag(curr_pos)])
    ax_vel.set_title(f'Velocity Field at t ≈ {field_times[frame]:.1f} s')
    return [quiv, scat_vel]

ani_vel = FuncAnimation(fig_vel, update_vel, frames=num_frames, init_func=init_vel, blit=False)
ani_vel.save('velocity_field_animation.mp4', writer='ffmpeg', fps=10)
plt.show()

# Animation 3: Velocity magnitude heatmap over time
fig_mag, ax_mag = plt.subplots(figsize=(12, 6))
ax_mag.plot(x_cyl, y_cyl, 'k-', linewidth=2)
ax_mag.set_xlim(x_min, x_max)
ax_mag.set_ylim(y_min, y_max)
ax_mag.set_xlabel('x')
ax_mag.set_ylabel('y')
ax_mag.set_title('Velocity Magnitude Heatmap Over Time')
ax_mag.grid(True)

# Initial pcolormesh
u0, v0, _ = velocity_fields[0]
mag0 = np.sqrt(u0**2 + v0**2)
pcm_mag = ax_mag.pcolormesh(X, Y, mag0, cmap='viridis', shading='auto', vmin=0, vmax=max_mag)
plt.colorbar(pcm_mag, ax=ax_mag, label='Velocity Magnitude')
scat_mag = ax_mag.scatter([], [], c='r', s=20, label='Vortices')

def init_mag():
    pcm_mag.set_array(np.zeros_like(mag0).ravel())
    scat_mag.set_offsets(np.empty((0, 2)))
    return [pcm_mag, scat_mag]

def update_mag(frame):
    u, v, _ = velocity_fields[frame]
    mag = np.sqrt(u**2 + v**2)
    pcm_mag.set_array(mag.ravel())
    curr_pos = positions_at_fields[frame]
    scat_mag.set_offsets(np.c_[np.real(curr_pos), np.imag(curr_pos)])
    ax_mag.set_title(f'Velocity Magnitude Heatmap at t ≈ {field_times[frame]:.1f} s')
    return [pcm_mag, scat_mag]

ani_mag = FuncAnimation(fig_mag, update_mag, frames=num_frames, init_func=init_mag, blit=False)
ani_mag.save('velocity_magnitude_animation.mp4', writer='ffmpeg', fps=10)
plt.show()

# Animation 4: Vorticity heatmap over time
fig_vort, ax_vort = plt.subplots(figsize=(12, 6))
ax_vort.plot(x_cyl, y_cyl, 'k-', linewidth=2)
ax_vort.set_xlim(x_min, x_max)
ax_vort.set_ylim(y_min, y_max)
ax_vort.set_xlabel('x')
ax_vort.set_ylabel('y')
ax_vort.set_title('Vorticity Heatmap Over Time')
ax_vort.grid(True)

# Initial pcolormesh
_, _, omega0 = velocity_fields[0]
pcm_vort = ax_vort.pcolormesh(X, Y, omega0, cmap='seismic', shading='auto', vmin=-max_omega, vmax=max_omega)
plt.colorbar(pcm_vort, ax=ax_vort, label='Vorticity')
scat_vort = ax_vort.scatter([], [], c='r', s=20, label='Vortices')

def init_vort():
    pcm_vort.set_array(np.zeros_like(omega0).ravel())
    scat_vort.set_offsets(np.empty((0, 2)))
    return [pcm_vort, scat_vort]

def update_vort(frame):
    _, _, omega = velocity_fields[frame]
    pcm_vort.set_array(omega.ravel())
    curr_pos = positions_at_fields[frame]
    scat_vort.set_offsets(np.c_[np.real(curr_pos), np.imag(curr_pos)])
    ax_vort.set_title(f'Vorticity Heatmap at t ≈ {field_times[frame]:.1f} s')
    return [pcm_vort, scat_vort]

ani_vort = FuncAnimation(fig_vort, update_vort, frames=num_frames, init_func=init_vort, blit=False)
ani_vort.save('vorticity_animation.mp4', writer='ffmpeg', fps=10)
plt.show()