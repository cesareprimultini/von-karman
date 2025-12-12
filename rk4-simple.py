"""Viscous vortex method for Von Kármán vortex streets shed from one or more cylinders."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from tqdm import tqdm
from scipy.interpolate import interp1d
from numba import njit, prange


# ============================================================================
# CONFIGURATION - All simulation parameters in one place
# ============================================================================

# --- Velocity configuration ---
velocity_mode = 'function'  # 'constant' or 'file'
U_inf_constant = 1.0  # Freestream velocity [m/s] (used if mode='constant')
velocity_file = 'velocity.xlsx'  # Path to Excel file (time in col A, velocity in col B)

# --- Geometry ---
D = 4.88  # Cylinder diameter [m]
spacing = 40.0  # Center-to-center spacing [m]

# Cylinder layout (easily configurable)
s = spacing / 2.0
cylinders = [
    {'x': -s, 'y': +s, 'D': D},
    {'x': +s, 'y': +s, 'D': D},
    {'x': -s, 'y': -s, 'D': D},
    {'x': +s, 'y': -s, 'D': D},
]

# --- Flow conditions ---
Re_initial = 1000000.0  # Reynolds number at t=0 (100=laminar, 6.7e6=transcritical offshore)
flow_angle_metocean = 270.0  # Flow direction [degrees, metocean convention]
rotation_angle = 30.0  # Cylinder cluster rotation [degrees]

# --- Turbulence modeling ---
enable_eddy_viscosity = (Re_initial > 1000)  # Enable turbulent diffusion for Re > 1000
enable_stochastic_shedding = (Re_initial > 1000)  # Enable random fluctuations for Re > 1000
enable_core_saturation = (Re_initial > 1000)  # Limit max core size for Re > 1000
sigma_max_factor = 0.5  # Maximum core size = sigma_max_factor * D_ref

# --- Time integration ---
dt = 0.01  # Time step [s]
total_time = 300.0  # Total simulation time [s]

# --- Vortex parameters ---
sigma_0_factor = 0.1  # Initial core size = sigma_0_factor * D_ref
theta_sep_deg = 80.0  # Separation angle [degrees] from rear stagnation point
x_removal = 250.0  # Remove vortices beyond this x-coordinate [m]

# --- Measurement probe ---
measure_point = (32.0, 5.0)  # (x, y) location [m]

# --- Visualization settings ---
plot_grid_size = 600  # Grid resolution for velocity field
plot_x_range = (-80.0, 150.0)  # (x_min, x_max) [m]
plot_y_range = (-50.0, 50.0)  # (y_min, y_max) [m]
plot_arrow_skip = 4  # Arrow decimation factor
save_wake_plot = True  # Save wake velocity field plot
save_velocity_history = True  # Save velocity time history plot
save_freestream_plot = True  # Save freestream velocity plot (if time-varying)
wake_plot_dpi = 600  # DPI for wake plot
history_plot_dpi = 150  # DPI for time history plots

# --- Derived parameters ---
D_ref = max(cyl['D'] for cyl in cylinders)  # Reference diameter [m]
flow_angle = 270.0 - flow_angle_metocean  # Internal flow angle [degrees]
theta_sep = np.radians(theta_sep_deg)  # Separation angle [radians]
sigma_0 = sigma_0_factor * D_ref  # Initial vortex core size [m]


# ============================================================================
# VELOCITY PROFILE SETUP
# ============================================================================

def load_velocity_from_file(filepath):
    """Load time-velocity data from Excel (column A = time, column B = velocity)"""
    import pandas as pd
    df = pd.read_excel(filepath)
    time_data = df.iloc[:, 0].values
    velocity_data = df.iloc[:, 1].values
    return interp1d(time_data, velocity_data, kind='linear', fill_value='extrapolate')


def velocity_ramp(t):
    """Example: Linear velocity ramp function"""
    return 1.0 + 0.01 * t


# Initialize velocity profile based on mode
if velocity_mode == 'file':
    velocity_interpolator = load_velocity_from_file(velocity_file)
    def get_velocity(t):
        """Returns freestream velocity at time t"""
        return float(velocity_interpolator(t))
elif velocity_mode == 'function':
    def get_velocity(t):
        """Returns freestream velocity at time t"""
        return velocity_ramp(t)
else:  # constant
    def get_velocity(t):
        """Returns freestream velocity at time t"""
        return U_inf_constant

# Compute fixed molecular viscosity from initial conditions
U_inf_initial = get_velocity(0.0)
nu = U_inf_initial * D_ref / Re_initial  # Kinematic viscosity [m²/s] (fixed, Re varies with U)


# ============================================================================
# PHYSICS FUNCTIONS
# ============================================================================

def compute_strouhal_number(Re):
    """Strouhal number from Reynolds (Roshko 1961 correlations extended for high Re)"""
    if Re < 47:
        return 0.0  # No vortex shedding
    elif Re < 300:
        return 0.212 * (1 - 21.2 / Re)  # Laminar regime
    elif Re < 2e5:
        return 0.20  # Subcritical regime
    elif Re < 3.5e6:
        log_Re = np.log10(Re)
        return 0.20 + 0.06 * (log_Re - 5.3)  # Critical/supercritical transition
    else:
        return 0.30  # Transcritical regime


def compute_circulation_magnitude(U_inf, D, St):
    """Target circulation Gamma = 2*pi*St*U*D (potential flow theory)"""
    return 2.0 * np.pi * St * U_inf * D


def compute_shedding_period(U_inf, D, St):
    """Shedding period T = D/(St*U) from Strouhal relation f = St*U/D"""
    if St <= 0:
        return np.inf
    return D / (St * U_inf)


def compute_eddy_viscosity_field(x, y, D, nu_molecular, U_inf, Re):
    """Effective viscosity (molecular + turbulent eddy contribution in wake)"""
    # Turbulent viscosity coefficient (calibrated from experiments)
    if Re < 1000:
        C_eddy = 0.01
    elif Re < 1e5:
        C_eddy = 0.02 + 0.01 * np.log10(Re / 1000)
    else:
        C_eddy = 0.10

    # Normalized coordinates
    x_norm = x / D
    y_norm = y / D

    # Wake profile: Gaussian transverse, exponential decay streamwise
    sigma_wake = 0.5 + 0.1 * x_norm  # Wake width grows linearly
    f_transverse = np.exp(-(y_norm**2) / (2 * sigma_wake**2))
    f_streamwise = np.exp(-x_norm / 20.0)  # Decays over ~20 diameters
    f_wake = f_transverse * f_streamwise

    # Total effective viscosity
    nu_turbulent = C_eddy * U_inf * D * f_wake
    return nu_molecular + nu_turbulent


def compute_formation_time(D, U_inf, Re):
    """Vortex formation time scale (formation number F* ~ 4.0)"""
    F_star = 4.0
    if Re > 1e5:
        F_star *= 0.8  # 20% faster formation in turbulent regime
    return F_star * (D / U_inf)


def circulation_growth_sigmoid(t_age, T_form):
    """Smooth circulation growth factor [0,1] during vortex formation"""
    if t_age < 0:
        return 0.0

    tau = t_age / T_form  # Normalized time
    k = 10.0  # Steepness

    # Sigmoid centered at tau = 0.5
    sigmoid_raw = 1.0 / (1.0 + np.exp(-k * (tau - 0.5)))

    # Normalize to [0, 1] range
    g_min = 1.0 / (1.0 + np.exp(k * 0.5))
    g_max = 1.0 / (1.0 + np.exp(-k * 0.5))
    g_normalized = (sigmoid_raw - g_min) / (g_max - g_min)

    return np.clip(g_normalized, 0.0, 1.0)


# ============================================================================
# VORTEX DYNAMICS
# ============================================================================

def apply_transforms(x, y, rot_deg, flow_deg):
    """Apply rotation and flow angle transforms to coordinates"""
    total_angle = np.radians(rot_deg + flow_deg)
    c, s = np.cos(total_angle), np.sin(total_angle)
    return c*x - s*y, s*x + c*y


def shed_vortex_with_turbulence(cyl, t, upper, Re, U_inf, D_ref, Gamma_mag, sigma_0,
                                  rotation_angle, flow_angle, theta_sep, enable_stochastic=True):
    """Create new vortex at separation point with turbulent fluctuations"""
    sign = 1 if upper else -1
    gamma_base = -Gamma_mag if upper else Gamma_mag

    # Add turbulent fluctuation to circulation
    if enable_stochastic and Re > 1000:
        fluctuation_level = 0.10 * np.log10(Re / 1000)
        fluctuation_level = min(fluctuation_level, 0.20)  # Cap at 20%
        gamma = gamma_base * (1 + np.random.normal(0, fluctuation_level))
    else:
        gamma = gamma_base

    # Compute separation angle with optional randomness
    if enable_stochastic and Re > 1000:
        theta_sep_perturbed = theta_sep + np.random.normal(0, np.radians(5))
    else:
        theta_sep_perturbed = theta_sep

    # Transform cylinder center to global coordinates
    x_cyl_global, y_cyl_global = apply_transforms(cyl['x'], cyl['y'], rotation_angle, flow_angle)

    # Shedding angle in global coordinates relative to flow direction
    flow_direction_rad = np.radians(flow_angle)
    rear_angle = flow_direction_rad + np.pi  # Rear stagnation point faces into flow
    theta = rear_angle - sign * theta_sep_perturbed

    # Shedding position in global coordinates
    a = cyl['D'] / 2
    x_shed = x_cyl_global + a * np.cos(theta)
    y_shed = y_cyl_global + a * np.sin(theta)

    return {
        'x': x_shed,
        'y': y_shed,
        'gamma_target': gamma,  # Maximum circulation (final value)
        'gamma': 0.0,  # Current circulation (starts at zero)
        'sigma': sigma_0,
        'birth_t': t,
        'U_inf_birth': U_inf,  # Velocity at birth time
        'T_form': compute_formation_time(D_ref, U_inf, Re)
    }


def remove_internal_vortices(cylinders, vortices, rotation_angle, flow_angle):
    """Filter vortices inside cylinders"""
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


# ============================================================================
# NUMBA KERNELS - DO NOT MODIFY (optimized for parallel performance)
# ============================================================================

@njit(parallel=True, fastmath=True)
def compute_velocities_numba(targets_x, targets_y, sources_x, sources_y, gammas, sigmas, U_inf_x, U_inf_y):
    """Numba kernel: velocity from vortices at target points"""
    M = len(targets_x)
    N = len(sources_x)
    U = np.empty(M, dtype=np.float64)
    V = np.empty(M, dtype=np.float64)

    for j in prange(M):
        u_total = U_inf_x
        v_total = U_inf_y

        for i in range(N):
            DX = targets_x[j] - sources_x[i]
            DY = targets_y[j] - sources_y[i]
            R_sq = DX * DX + DY * DY

            if R_sq > 1e-12:  # Avoid singularity
                R = np.sqrt(R_sq)
                sigma_sq = sigmas[i] * sigmas[i]
                EXP = np.exp(-R_sq / sigma_sq)
                V_theta = (gammas[i] / (2.0 * np.pi * R)) * (1.0 - EXP)

                u_total += -V_theta * DY / R
                v_total += V_theta * DX / R

        U[j] = u_total
        V[j] = v_total

    return U, V


@njit(parallel=True, fastmath=True)
def compute_velocity_field_numba(X_flat, Y_flat, vortex_x, vortex_y, vortex_gamma, vortex_sigma, U_inf_x, U_inf_y):
    """Numba kernel: velocity field over meshgrid"""
    M = len(X_flat)
    N = len(vortex_x)
    U = np.empty(M, dtype=np.float64)
    V = np.empty(M, dtype=np.float64)

    for j in prange(M):
        u_total = U_inf_x
        v_total = U_inf_y

        for i in range(N):
            DX = X_flat[j] - vortex_x[i]
            DY = Y_flat[j] - vortex_y[i]
            R_sq = DX * DX + DY * DY

            if R_sq > 1e-12:
                R = np.sqrt(R_sq)
                sigma_sq = vortex_sigma[i] * vortex_sigma[i]
                EXP = np.exp(-R_sq / sigma_sq)
                V_theta = (vortex_gamma[i] / (2.0 * np.pi * R)) * (1.0 - EXP)

                u_total += -V_theta * DY / R
                v_total += V_theta * DX / R

        U[j] = u_total
        V[j] = v_total

    return U, V


# ============================================================================
# VELOCITY FIELD COMPUTATION
# ============================================================================

def compute_velocities_at_points(targets_x, targets_y, sources_x, sources_y, gammas, sigmas, U_inf, flow_angle):
    """Calculate velocity at arbitrary points (no boundary enforcement, transparent cylinders)"""
    M = len(targets_x)
    if M == 0:
        return np.array([]), np.array([])

    # Compute freestream components
    flow_angle_rad = np.radians(flow_angle)
    U_inf_x = U_inf * np.cos(flow_angle_rad)
    U_inf_y = U_inf * np.sin(flow_angle_rad)

    # Call optimized numba kernel
    return compute_velocities_numba(targets_x, targets_y, sources_x, sources_y, gammas, sigmas, U_inf_x, U_inf_y)


def compute_velocity_field(X, Y, vortices, U_inf, flow_angle):
    """Compute velocity field over meshgrid"""
    # Extract vortex properties
    vortex_x = np.array([v['x'] for v in vortices])
    vortex_y = np.array([v['y'] for v in vortices])
    vortex_gamma = np.array([v['gamma'] for v in vortices])
    vortex_sigma = np.array([v['sigma'] for v in vortices])

    # Compute freestream components
    flow_angle_rad = np.radians(flow_angle)
    U_inf_x = U_inf * np.cos(flow_angle_rad)
    U_inf_y = U_inf * np.sin(flow_angle_rad)

    # Flatten meshgrid for numba processing
    X_flat = X.ravel()
    Y_flat = Y.ravel()

    # Call numba kernel
    U_flat, V_flat = compute_velocity_field_numba(
        X_flat, Y_flat, vortex_x, vortex_y, vortex_gamma, vortex_sigma, U_inf_x, U_inf_y
    )

    # Reshape back to meshgrid
    U = U_flat.reshape(X.shape)
    V = V_flat.reshape(X.shape)

    return U, V


# ============================================================================
# MAIN SIMULATION
# ============================================================================

# Initialize simulation
all_vortices = []
time = 0.0
shed_counters = [0 for _ in cylinders]

# Initialize shedding schedule with initial velocity
U_inf_0 = get_velocity(0.0)
Re_0 = U_inf_0 * D_ref / nu
St_0 = compute_strouhal_number(Re_0)
shed_period_0 = compute_shedding_period(U_inf_0, D_ref, St_0)
next_shed_times = [shed_period_0 for _ in cylinders]

# Data storage (easy to extend with new variables)
history_data = {
    'time': [],
    'ux': [],
    'uy': [],
    'vmag': []
}
shed_times = []  # Track shedding events

# Print initial conditions
print("Starting vortex method simulation...")
print(f"  U_inf(0) = {U_inf_0:.3f} m/s")
print(f"  Re(0) = {Re_0:.2e}")
print(f"  St(0) = {St_0:.3f}")
print(f"  nu (fixed) = {nu:.6e} m²/s")
print(f"  Initial shedding period = {shed_period_0:.3f} s\n")

num_steps = int(total_time / dt)

# Time integration loop
for step in tqdm(range(num_steps), desc="Simulating", mininterval=0.5, unit="step"):

    # 1. VORTEX SHEDDING
    for c_idx, cyl in enumerate(cylinders):
        if time >= next_shed_times[c_idx]:
            # Evaluate instantaneous velocity and conditions
            U_inf_now = get_velocity(time)
            Re_now = U_inf_now * D_ref / nu
            St_now = compute_strouhal_number(Re_now)
            Gamma_mag_now = compute_circulation_magnitude(U_inf_now, D_ref, St_now)

            # Shed vortex with current conditions
            upper = (shed_counters[c_idx] % 2 == 0)
            new_vortex = shed_vortex_with_turbulence(
                cyl, time, upper, Re_now, U_inf_now, D_ref, Gamma_mag_now, sigma_0,
                rotation_angle, flow_angle, theta_sep, enable_stochastic=enable_stochastic_shedding
            )
            all_vortices.append(new_vortex)
            shed_counters[c_idx] += 1

            # Schedule next shedding with current velocity
            shed_period_now = compute_shedding_period(U_inf_now, D_ref, St_now)
            next_shed_times[c_idx] += shed_period_now
            shed_times.append(time)

    N = len(all_vortices)
    if N == 0:
        time += dt
        continue

    # 2. CORE SIZE UPDATE (viscous diffusion)
    birth_times = np.array([v['birth_t'] for v in all_vortices])
    ages = time - birth_times

    if enable_eddy_viscosity:
        # Spatially-varying effective viscosity
        U_inf_now = get_velocity(time)
        Re_now = U_inf_now * D_ref / nu

        for i in range(N):
            x_vortex = all_vortices[i]['x']
            y_vortex = all_vortices[i]['y']
            nu_eff = compute_eddy_viscosity_field(x_vortex, y_vortex, D_ref, nu, U_inf_now, Re_now)
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

    # 3. CIRCULATION GROWTH (gradual formation)
    for i in range(N):
        t_age = time - all_vortices[i]['birth_t']
        growth_factor = circulation_growth_sigmoid(t_age, all_vortices[i]['T_form'])
        all_vortices[i]['gamma'] = all_vortices[i]['gamma_target'] * growth_factor

    # 4. RK4 ADVECTION (time-accurate velocity evaluation)
    sources_x = np.array([v['x'] for v in all_vortices])
    sources_y = np.array([v['y'] for v in all_vortices])
    gammas = np.array([v['gamma'] for v in all_vortices])
    sigmas = np.array([v['sigma'] for v in all_vortices])

    # Stage 1: t
    U_inf_t = get_velocity(time)
    U1, V1 = compute_velocities_at_points(sources_x, sources_y, sources_x, sources_y, gammas, sigmas, U_inf_t, flow_angle)
    k1x = dt * U1
    k1y = dt * V1

    # Stage 2: t + 0.5*dt
    U_inf_mid = get_velocity(time + 0.5*dt)
    mid_x = sources_x + 0.5 * k1x
    mid_y = sources_y + 0.5 * k1y
    U2, V2 = compute_velocities_at_points(mid_x, mid_y, sources_x, sources_y, gammas, sigmas, U_inf_mid, flow_angle)
    k2x = dt * U2
    k2y = dt * V2

    # Stage 3: t + 0.5*dt
    mid_x = sources_x + 0.5 * k2x
    mid_y = sources_y + 0.5 * k2y
    U3, V3 = compute_velocities_at_points(mid_x, mid_y, sources_x, sources_y, gammas, sigmas, U_inf_mid, flow_angle)
    k3x = dt * U3
    k3y = dt * V3

    # Stage 4: t + dt
    U_inf_end = get_velocity(time + dt)
    end_x = sources_x + k3x
    end_y = sources_y + k3y
    U4, V4 = compute_velocities_at_points(end_x, end_y, sources_x, sources_y, gammas, sigmas, U_inf_end, flow_angle)
    k4x = dt * U4
    k4y = dt * V4

    # Update positions
    new_x = sources_x + (k1x + 2*k2x + 2*k3x + k4x) / 6
    new_y = sources_y + (k1y + 2*k2y + 2*k3y + k4y) / 6
    for i in range(N):
        all_vortices[i]['x'] = new_x[i]
        all_vortices[i]['y'] = new_y[i]

    # 5. VORTEX REMOVAL
    all_vortices = [v for v in all_vortices if v['x'] <= x_removal]

    # 6. MEASUREMENTS
    U_inf_now = get_velocity(time)
    ux, uy = compute_velocities_at_points(
        np.array([measure_point[0]]),
        np.array([measure_point[1]]),
        sources_x, sources_y, gammas, sigmas, U_inf_now, flow_angle
    )
    history_data['time'].append(time)
    history_data['ux'].append(ux[0] if len(ux) > 0 else U_inf_now * np.cos(np.radians(flow_angle)))
    history_data['uy'].append(uy[0] if len(uy) > 0 else U_inf_now * np.sin(np.radians(flow_angle)))
    history_data['vmag'].append(np.sqrt(history_data['ux'][-1]**2 + history_data['uy'][-1]**2))

    time += dt

print(f"Simulation complete. Total vortices in domain: {len(all_vortices)}")


# ============================================================================
# VISUALIZATION
# ============================================================================

# Compute velocity field for visualization
if save_wake_plot:
    print("Computing velocity field...")
    x = np.linspace(plot_x_range[0], plot_x_range[1], plot_grid_size)
    y = np.linspace(plot_y_range[0], plot_y_range[1], plot_grid_size)
    X, Y = np.meshgrid(x, y)

    U_inf_final = get_velocity(time)
    U, V = compute_velocity_field(X, Y, all_vortices, U_inf_final, flow_angle)
    vel_mag = np.sqrt(U**2 + V**2)

    # Plot wake velocity field
    fig, ax = plt.subplots(figsize=(16, 8), dpi=300)

    contour = ax.contourf(X, Y, vel_mag, levels=100, cmap='viridis')
    plt.colorbar(contour, ax=ax, label='Velocity Magnitude [m/s]')

    # Velocity arrows
    skip = plot_arrow_skip
    scale_reference = np.percentile(vel_mag, 95)
    scale_dynamic = scale_reference * 80
    ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
              U[::skip, ::skip], V[::skip, ::skip],
              color='white', alpha=0.7, scale=scale_dynamic, width=0.0005)

    # Draw cylinders
    for cyl in cylinders:
        x_cyl, y_cyl = apply_transforms(cyl['x'], cyl['y'], rotation_angle, flow_angle)
        circle = Circle((x_cyl, y_cyl), cyl['D']/2,
                        color='gray', fill=True, zorder=10,
                        edgecolor='black', linewidth=2)
        ax.add_patch(circle)

    # Draw vortices
    for v in all_vortices:
        color = 'ro' if v['gamma'] > 0 else 'bo'
        ax.plot(v['x'], v['y'], color, markersize=2, alpha=0.4)

    # Draw measurement probe
    ax.plot(measure_point[0], measure_point[1], 'r*', markersize=12,
            markeredgewidth=1.5, markeredgecolor='black', label='CPS', zorder=15)

    ax.set_xlim(plot_x_range)
    ax.set_ylim(plot_y_range)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    Re_final = U_inf_final * D_ref / nu
    ax.set_title(f'Vortex Wake (D={D_ref}m, Re={Re_final:.2e}, t={total_time}s)')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig('cylinder_wake.png', dpi=wake_plot_dpi, bbox_inches='tight')
    print("Saved cylinder_wake.png")
    plt.close()


# Velocity time history plot
if save_velocity_history:
    print("Generating velocity measurement plot...")
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(history_data['time'], history_data['ux'], label='$u_x$ (streamwise)',
            linewidth=1.5, alpha=0.8)
    ax.plot(history_data['time'], history_data['uy'], label='$u_y$ (transverse)',
            linewidth=1.5, alpha=0.8)
    ax.plot(history_data['time'], history_data['vmag'], label='$|v|$ (magnitude)',
            linewidth=2, color='black')

    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Velocity [m/s]', fontsize=12)
    ax.set_title(f'Velocity Time History at ({measure_point[0]:.1f}, {measure_point[1]:.1f}) m',
                 fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11)
    ax.set_xlim([0, total_time])

    # Mark shedding events
    for t_shed in shed_times:
        ax.axvline(t_shed, color='red', alpha=0.15, linewidth=0.5, linestyle='-')

    plt.tight_layout()
    plt.savefig('velocity_measurement.png', dpi=history_plot_dpi, bbox_inches='tight')
    print("Saved velocity_measurement.png")
    plt.close()


# Freestream velocity plot (if time-varying)
if velocity_mode != 'constant' and save_freestream_plot:
    print("Generating freestream velocity plot...")
    fig, ax = plt.subplots(figsize=(10, 4))

    t_plot = np.array(history_data['time'])
    U_plot = np.array([get_velocity(t) for t in t_plot])
    Re_plot = U_plot * D_ref / nu

    ax.plot(t_plot, U_plot, 'b-', linewidth=2, label='$U_\\infty(t)$')
    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Freestream Velocity [m/s]', fontsize=12, color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.set_title('Time-Varying Freestream Velocity', fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, total_time])

    # Secondary axis for Reynolds number
    ax2 = ax.twinx()
    ax2.plot(t_plot, Re_plot, 'r--', linewidth=2, label='Re(t)')
    ax2.set_ylabel('Reynolds Number', fontsize=12, color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    ax2.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

    # Combined legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=11)

    plt.tight_layout()
    plt.savefig('freestream_velocity_timeseries.png', dpi=history_plot_dpi, bbox_inches='tight')
    print("Saved freestream_velocity_timeseries.png")
    plt.close()


print(f"\nDone! S/D = {spacing/D_ref:.1f}")
print("This vortex wake field can be combined with:")
print("  1. Potential flow (local cylinder effects)")
print("  2. Global flow (terrain/hill effects)")
print("  for total velocity amplification at probe locations")
