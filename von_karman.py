"""
von Kármán Vortex Street Simulator
Pure functional implementation of algebraic point vortex model
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from dataclasses import dataclass

# ============================================================================
# VORTEX DATA STRUCTURE (Dynamic Tracking)
# ============================================================================

@dataclass
class Vortex:
    """
    Represents a single vortex in the dynamic vortex street.

    Attributes:
        x: Streamwise position [D]
        y: Transverse position [D]
        gamma: Circulation strength (signed) [m²/s]
        t_birth: Time when vortex was shed [T = D/U∞]
    """
    x: float
    y: float
    gamma: float
    t_birth: float

    def sigma(self, t, params):
        """
        Compute time-dependent core radius via Lamb-Oseen formula.

        Physics: Vortex cores grow due to viscous diffusion
        Formula: σ²(t) = σ₀² + 4ν(t - t_birth)

        Args:
            t: Current time [T = D/U∞]
            params: Parameter dictionary

        Returns:
            Core radius σ(t) [D]
        """
        age = t - self.t_birth
        nu = params['U_inf'] * params['D'] / params['Re']
        sigma_sq = params['sigma_0']**2 + 4 * nu * age
        return np.sqrt(sigma_sq)

    def __repr__(self):
        return f"Vortex(x={self.x:.2f}, y={self.y:.2f}, Γ={self.gamma:+.3f})"

# ============================================================================
# PARAMETER CONFIGURATION
# ============================================================================

# Physics: von Kármán stability condition h/l ≈ 0.281
# Reynolds regimes: Re < 47 (no shedding), 47 < Re < 190 (stable), Re > 190 (turbulent)
params = {
    # Vortex strength and geometry
    'Gamma': 1.5,          # Circulation strength [m²/s or normalized]
    'h': 0.281,            # Initial transverse spacing between vortex rows [D]
    'l': 1.0,              # Streamwise period (wavelength) [D]
    'a': 0.5,              # Stagger parameter: 0.5 means a = l/2 (full stagger)
    'Re': 60.0,           # Reynolds number (stable laminar street)
    'D': 1.0,              # Cylinder diameter [m or normalized to 1]
    'U_inf': 1.0,          # Freestream velocity [m/s or normalized]

    # Viscous vortex blob model parameters (NEW!)
    'sigma_0': 0.1,        # Initial vortex core radius [D] - typical: 0.05-0.15
    'spreading_coef': 0.10, # Wake spreading coefficient α - typical: 0.01-0.05 (Re-dependent)
                            # TRY: 0.08-0.10 for MORE visible divergence
                            #      0.03-0.05 for subtle/realistic divergence
                            #      0.01-0.02 for minimal divergence (high Re)

    # Grid and domain
    'grid_nx': 250,        # Grid resolution in x (streamwise)
    'grid_ny': 120,        # Grid resolution in y (transverse)
    'x_min': -1.0,         # Domain minimum x [D]
    'x_max': 20.0,         # Domain maximum x [D]
    'y_min': -2.5,         # Domain minimum y [D]
    'y_max': 2.5,          # Domain maximum y [D]
    'n_vortices': 15,      # Number of vortex pairs (±n, total 2n)

    # Time-varying dynamics
    'probe_x': 4.0,        # Probe location for time series x [D]
    'probe_y': 0.3,        # Probe location for time series y [D]
    'n_timesteps': 60,     # Number of time steps for temporal evolution (DEPRECATED)

    # Decay model
    'decay_alpha': 0.08,   # Decay rate coefficient (dimensionless)
    'decay_x0': 1.5,       # Decay onset position [D]

    # ========================================================================
    # DYNAMIC VORTEX TRACKING PARAMETERS (NEW!)
    # ========================================================================
    'St': 0.2,              # Strouhal number (typical for Re~100)
    'C_gamma': 2 * np.pi,   # Circulation coefficient (TUNABLE: 4.0-8.0)
    'theta_sep': 1.4,       # Separation angle [rad] (~80 deg)
    't_max': 20.0,          # Total simulation time [T = D/U_inf] (reduced for testing)
    'dt': 0.02,             # Integration timestep [T] (was 0.05)
    'epsilon': 0.01,        # Self-induction cutoff radius [D]
    'use_image_vortices': True,  # Enable image vortex boundary conditions

    # Snapshot configuration
    'n_snapshots': 50,      # Number of velocity field snapshots to save
    'snapshots_folder': 'out',  # Folder for snapshot images
}

# ============================================================================
# GEOMETRY FUNCTIONS
# ============================================================================

def get_vortex_positions(params):
    """
    Generate staggered vortex street positions with DIVERGING rows (viscous model).

    Physics: Viscous vortex blob model with wake spreading
    - Upper row: +Γ at (k·l, +h(x)/2) - DIVERGES downstream
    - Lower row: -Γ at (a + k·l, -h(x)/2) - DIVERGES downstream
    - Each vortex has core radius σ(x) that grows with x

    Key differences from inviscid model:
    - Transverse spacing h(x) INCREASES downstream (not constant)
    - Core radii σ(x) GROW downstream via Lamb-Oseen formula
    - Rows DIVERGE (not parallel)

    Stagger a = l/2 gives classical von Kármán configuration
    """
    l = params['l']
    a = params['a'] * l  # Convert normalized to absolute
    n = params['n_vortices']

    # Upper row: positive circulation with DIVERGING positions
    x_pos = np.array([k * l for k in range(-n, n + 1)])
    # Y-position varies with x due to wake spreading
    y_pos = np.array([compute_transverse_spacing(x, params) / 2 for x in x_pos])
    # Core radius varies with x due to viscous diffusion
    sigma_pos = np.array([compute_vortex_core_radius(x, params) for x in x_pos])

    # Lower row: negative circulation with DIVERGING positions
    x_neg = np.array([a + k * l for k in range(-n, n + 1)])
    # Y-position varies with x (symmetric about centerline)
    y_neg = np.array([-compute_transverse_spacing(x, params) / 2 for x in x_neg])
    # Core radius varies with x
    sigma_neg = np.array([compute_vortex_core_radius(x, params) for x in x_neg])

    return {
        'x_pos': x_pos, 'y_pos': y_pos, 'sigma_pos': sigma_pos,
        'x_neg': x_neg, 'y_neg': y_neg, 'sigma_neg': sigma_neg
    }


def compute_cone_angle(params):
    """
    Calculate cone opening angle from vortex spacing geometry.

    Physics: Cone boundary connects successive vortices
    Geometry: tan(θ) = h / (l/2) for staggered rows
    """
    h = params['h']
    l = params['l']

    # Half-angle: opposite = h, adjacent = l/2
    theta_half_rad = np.arctan(h / (l / 2))
    theta_half_deg = np.degrees(theta_half_rad)
    theta_full_deg = 2 * theta_half_deg

    return {
        'theta_half_rad': theta_half_rad,
        'theta_half_deg': theta_half_deg,
        'theta_full_deg': theta_full_deg
    }


def get_cone_boundaries(params, x_range):
    """
    Return upper and lower cone boundary coordinates.

    Physics: Cone emanates from cylinder with half-angle θ
    Boundaries: y = ±x·tan(θ)
    """
    theta = compute_cone_angle(params)['theta_half_rad']
    y_upper = x_range * np.tan(theta)
    y_lower = -x_range * np.tan(theta)
    return y_upper, y_lower


# ============================================================================
# IMAGE VORTEX METHOD (Cylinder Boundary Conditions)
# ============================================================================

def get_image_vortex_position(vortex, params):
    """
    Compute image vortex position for cylinder boundary condition.

    Physics: Image vortices ensure no flow through cylinder surface (no-penetration BC)
    For vortex at (x,y), image is at (x_img, y_img) where:
        z_img = a²/z̄  (complex conjugate)
    In Cartesian:
        x_img = a²x/(x²+y²)
        y_img = -a²y/(x²+y²)

    The image has OPPOSITE circulation to satisfy the boundary condition.

    Args:
        vortex: Vortex object
        params: Parameter dictionary

    Returns:
        tuple (x_img, y_img) or None if vortex inside cylinder
    """
    a = params['D'] / 2.0  # Cylinder radius
    r_sq = vortex.x**2 + vortex.y**2

    # Avoid singularity if vortex is inside or very close to cylinder
    if r_sq < (a * 1.1)**2:
        return None  # No valid image

    x_img = a**2 * vortex.x / r_sq
    y_img = -a**2 * vortex.y / r_sq  # Negative for complex conjugate

    return x_img, y_img


# ============================================================================
# VORTEX SHEDDING LOGIC
# ============================================================================

def shed_vortex(vortices, t, params, shedding_state):
    """
    Shed vortex at regular intervals based on Strouhal number.

    Shedding alternates between upper (negative Γ, clockwise) and
    lower (positive Γ, counter-clockwise) with positions based on
    separation angle on cylinder surface.

    Physics:
    - Shedding frequency: f = St·U∞/D
    - Period: T_shed = D/(St·U∞)
    - Circulation: Γ = C_gamma · St · U∞ · D (tunable via C_gamma)
    - Position: On cylinder at separation angle θ_sep

    Args:
        vortices: List of Vortex objects (modified in-place)
        t: Current time [T = D/U∞]
        params: Parameter dictionary
        shedding_state: Dictionary with 'count' key tracking number shed

    Sign convention (follows blevins_vortex.py):
    - Upper vortex: negative Γ (clockwise rotation)
    - Lower vortex: positive Γ (counter-clockwise rotation)
    """
    St = params['St']
    D = params['D']
    U_inf = params['U_inf']
    C_gamma = params['C_gamma']
    a = D / 2.0  # Cylinder radius
    theta_sep = params['theta_sep']

    # Shedding period (time between consecutive vortices)
    T_shed = D / (St * U_inf)

    # How many vortices should have been shed by now?
    expected_count = int(t / T_shed)

    # Shed vortices if we're behind schedule
    while shedding_state['count'] < expected_count:
        shed_idx = shedding_state['count']

        # Alternate between upper and lower
        is_upper = (shed_idx % 2 == 0)

        # Position on cylinder surface at separation angle
        # Upper: π - θ_sep (left side, upper)
        # Lower: π + θ_sep (left side, lower)
        theta = np.pi - theta_sep if is_upper else np.pi + theta_sep

        # Shed slightly outside cylinder surface (at radius a + δ)
        delta_r = params['sigma_0']  # Offset by initial core size
        x_shed = (a + delta_r) * np.cos(theta)
        y_shed = (a + delta_r) * np.sin(theta)

        # Circulation magnitude (TUNABLE via C_gamma parameter)
        Gamma_shed = C_gamma * St * U_inf * D

        # Sign convention: upper=clockwise (negative), lower=CCW (positive)
        gamma_signed = -Gamma_shed if is_upper else Gamma_shed

        # Create new vortex
        new_vortex = Vortex(
            x=x_shed,
            y=y_shed,
            gamma=gamma_signed,
            t_birth=shed_idx * T_shed
        )

        vortices.append(new_vortex)
        shedding_state['count'] += 1

        # Debug output for first few vortices
        if shed_idx < 5:
            print(f"  Shed vortex {shed_idx}: {'UPPER' if is_upper else 'LOWER'} "
                  f"at ({x_shed:.3f}, {y_shed:.3f}), Gamma={gamma_signed:+.3f}")


# ============================================================================
# VORTEX ADVECTION (RK4 Integration)
# ============================================================================

def compute_vortex_velocity(vortex, vortices, t, params, exclude_idx=None):
    """
    Compute velocity at a vortex position due to all other vortices + images.

    Physics: dx/dt = u(x, t) where u is the velocity field
    This implements the Biot-Savart-like induction with Lamb-Oseen cores.

    Args:
        vortex: Vortex object whose velocity we're computing
        vortices: List of all vortices
        t: Current time [T = D/U∞]
        params: Parameter dictionary
        exclude_idx: Index to exclude (avoid self-induction)

    Returns:
        np.array([u, v]): Velocity components
    """
    epsilon = params['epsilon']
    D = params['D']
    U_inf = params['U_inf']
    a = D / 2.0

    u = U_inf  # Freestream contribution
    v = 0.0

    for i, other in enumerate(vortices):
        if exclude_idx is not None and i == exclude_idx:
            continue  # Skip self-induction

        # Direct vortex contribution
        dx = vortex.x - other.x
        dy = vortex.y - other.y
        r_sq = dx**2 + dy**2

        # Self-induction cutoff (avoid singularity when vortices are close)
        if r_sq < (epsilon * D)**2:
            continue

        r = np.sqrt(r_sq)
        sigma = other.sigma(t, params)
        exp_factor = np.exp(-r_sq / sigma**2)
        vel_mag = (other.gamma / (2 * np.pi * r)) * (1 - exp_factor)

        u += vel_mag * (-dy / r)
        v += vel_mag * (dx / r)

        # Image vortex contribution (for cylinder boundary condition)
        if params.get('use_image_vortices', False):
            image_pos = get_image_vortex_position(other, params)
            if image_pos is not None:
                x_img, y_img = image_pos
                dx_img = vortex.x - x_img
                dy_img = vortex.y - y_img
                r_sq_img = dx_img**2 + dy_img**2

                if r_sq_img >= (epsilon * D)**2:
                    r_img = np.sqrt(r_sq_img)
                    # Image has OPPOSITE circulation
                    vel_mag_img = (-other.gamma / (2 * np.pi * r_img)) * (
                        1 - np.exp(-r_sq_img / sigma**2)
                    )
                    u += vel_mag_img * (-dy_img / r_img)
                    v += vel_mag_img * (dx_img / r_img)

    return np.array([u, v])


def advect_vortices_rk4(vortices, t, dt, params):
    """
    4th-order Runge-Kutta integration for vortex positions.

    Integrates the ODE: dx/dt = u(x, t) for all vortices

    Algorithm (classical RK4):
        k1 = f(t, x0)
        k2 = f(t + dt/2, x0 + k1*dt/2)
        k3 = f(t + dt/2, x0 + k2*dt/2)
        k4 = f(t + dt, x0 + k3*dt)
        x_new = x0 + (dt/6)(k1 + 2k2 + 2k3 + k4)

    Stability: Requires dt < σ₀/U_max (typically dt ~ 0.02 for σ₀=0.1, U=1)

    Args:
        vortices: List of Vortex objects (modified in-place)
        t: Current time [T = D/U∞]
        dt: Timestep [T]
        params: Parameter dictionary
    """
    n = len(vortices)
    if n == 0:
        return

    # Store initial positions
    x0 = np.array([[v.x, v.y] for v in vortices])

    # Stage 1: k1 = f(t, x0)
    k1 = np.array([
        compute_vortex_velocity(vortices[i], vortices, t, params, exclude_idx=i)
        for i in range(n)
    ])

    # Stage 2: k2 = f(t + dt/2, x0 + k1*dt/2)
    for i in range(n):
        vortices[i].x, vortices[i].y = x0[i] + 0.5 * dt * k1[i]
    k2 = np.array([
        compute_vortex_velocity(vortices[i], vortices, t + 0.5*dt, params, exclude_idx=i)
        for i in range(n)
    ])

    # Stage 3: k3 = f(t + dt/2, x0 + k2*dt/2)
    for i in range(n):
        vortices[i].x, vortices[i].y = x0[i] + 0.5 * dt * k2[i]
    k3 = np.array([
        compute_vortex_velocity(vortices[i], vortices, t + 0.5*dt, params, exclude_idx=i)
        for i in range(n)
    ])

    # Stage 4: k4 = f(t + dt, x0 + k3*dt)
    for i in range(n):
        vortices[i].x, vortices[i].y = x0[i] + dt * k3[i]
    k4 = np.array([
        compute_vortex_velocity(vortices[i], vortices, t + dt, params, exclude_idx=i)
        for i in range(n)
    ])

    # Final weighted update: x_new = x0 + (dt/6)(k1 + 2k2 + 2k3 + k4)
    for i in range(n):
        pos_new = x0[i] + (dt / 6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
        vortices[i].x, vortices[i].y = pos_new


# ============================================================================
# VISCOUS VORTEX BLOB MODEL: LAMB-OSEEN AND WAKE SPREADING
# ============================================================================

def compute_vortex_core_radius(x, params):
    """
    Lamb-Oseen vortex core growth model (exact Navier-Stokes solution).

    Physics: Viscous diffusion causes vortex cores to grow as
    σ²(x) = σ₀² + 4νx/U∞ = σ₀² + 4xD/Re

    This is derived from the heat equation for vorticity diffusion.

    Args:
        x: Downstream distance (can be scalar or array) [D]
        params: Parameter dictionary

    Returns:
        Core radius σ(x) [D]
    """
    sigma_0 = params['sigma_0']  # Initial core radius
    Re = params['Re']
    D = params['D']

    # Handle negative x (upstream of cylinder) - use initial core size
    x_positive = np.maximum(x, 0)

    # Lamb-Oseen viscous spreading formula
    sigma_squared = sigma_0**2 + 4 * x_positive * D / Re

    return np.sqrt(sigma_squared)


def compute_transverse_spacing(x, params):
    """
    Wake width spreading model - empirical linear growth.

    Physics: Viscous diffusion and turbulent mixing cause wake to spread
    h(x) = h₀(1 + α·x/D)

    Empirical evidence shows linear growth in near wake, confirmed by
    experiments and CFD simulations.

    Args:
        x: Downstream distance (can be scalar or array) [D]
        params: Parameter dictionary

    Returns:
        Transverse spacing h(x) [D]
    """
    h_0 = params['h']  # Initial spacing (von Kármán stability: h/l ≈ 0.281)
    alpha = params['spreading_coef']  # Empirical spreading coefficient
    D = params['D']

    # Handle negative x (upstream) - use initial spacing
    x_positive = np.maximum(x, 0)

    # Linear wake spreading model
    h_x = h_0 * (1 + alpha * x_positive / D)

    return h_x


def gaussian_velocity_contribution(X, Y, x_v, y_v, sigma, Gamma):
    """
    Compute velocity contribution from a single Gaussian vortex blob.

    Physics: Lamb-Oseen vortex velocity field
    V_θ(r) = (Γ/2πr)[1 - exp(-(r/σ)²)]

    Unlike point vortices (singular at r=0), Gaussian blobs have:
    - Smooth velocity field (no singularities)
    - V → 0 as r → 0 (physically correct)
    - V → Γ/(2πr) as r → ∞ (matches point vortex far field)

    Args:
        X, Y: Meshgrid coordinates
        x_v, y_v: Vortex center position
        sigma: Vortex core radius
        Gamma: Circulation strength

    Returns:
        U, V: Velocity components
    """
    dx = X - x_v
    dy = Y - y_v
    r_squared = dx**2 + dy**2
    r = np.sqrt(r_squared + 1e-10)  # Small epsilon to avoid division by zero at center

    # Lamb-Oseen velocity magnitude with Gaussian cutoff
    exp_factor = np.exp(-r_squared / sigma**2)
    velocity_magnitude = (Gamma / (2 * np.pi * r)) * (1 - exp_factor)

    # Convert to Cartesian components (perpendicular to radius vector)
    # For 2D vortex: v = V_θ(r) * θ_hat = V_θ(r) * (-sin(θ), cos(θ))
    # where sin(θ) = dy/r, cos(θ) = dx/r
    U = velocity_magnitude * (-dy / r)
    V = velocity_magnitude * (dx / r)

    return U, V


def gaussian_velocity_contribution_dynamic(X, Y, vortex, t, params):
    """
    Compute velocity from single vortex with TIME-DEPENDENT core radius.

    This is the dynamic version where σ depends on vortex age, not position.

    Args:
        X, Y: Meshgrid coordinates or scalar positions
        vortex: Vortex object
        t: Current time [T = D/U∞]
        params: Parameter dictionary

    Returns:
        U, V: Velocity components
    """
    # Time-dependent core radius (Lamb-Oseen growth)
    sigma = vortex.sigma(t, params)

    # Use the original function with time-dependent sigma
    return gaussian_velocity_contribution(X, Y, vortex.x, vortex.y, sigma, vortex.gamma)


# ============================================================================
# PHYSICS: VELOCITY AND VORTICITY
# ============================================================================

def compute_velocity_field(X, Y, params):
    """
    Compute 2D velocity field using Gaussian vortex blobs (Lamb-Oseen model).

    Physics: Viscous vortex blob method replaces singular point vortices with
    Gaussian cores that grow due to viscous diffusion.

    V_θ(r, σ) = (Γ/2πr)[1 - exp(-(r/σ)²)]

    Key improvements over point vortex model:
    - No singularities (smooth velocity field)
    - Core size σ(x) varies with downstream distance
    - Physically realistic (exact Navier-Stokes solution)
    - Rows diverge downstream (matches experiments)
    """
    Gamma = params['Gamma']
    vortex_pos = get_vortex_positions(params)

    # Initialize velocity components
    U = np.zeros_like(X)
    V = np.zeros_like(X)

    # Contribution from positive circulation vortices (Gaussian blobs)
    for x_v, y_v, sigma in zip(vortex_pos['x_pos'], vortex_pos['y_pos'], vortex_pos['sigma_pos']):
        U_contrib, V_contrib = gaussian_velocity_contribution(X, Y, x_v, y_v, sigma, Gamma)
        U += U_contrib
        V += V_contrib

    # Contribution from negative circulation vortices (Gaussian blobs)
    for x_v, y_v, sigma in zip(vortex_pos['x_neg'], vortex_pos['y_neg'], vortex_pos['sigma_neg']):
        U_contrib, V_contrib = gaussian_velocity_contribution(X, Y, x_v, y_v, sigma, -Gamma)
        U += U_contrib
        V += V_contrib

    # Add freestream velocity
    U += params['U_inf']

    return U, V


def compute_vorticity_field(X, Y, params):
    """
    Compute vorticity field from velocity using finite differences.

    Physics: ω = ∂v/∂x - ∂u/∂y (curl of velocity in 2D)
    Vorticity concentrates at vortex cores

    Numerical: Central differences for 2nd order accuracy
    """
    U, V = compute_velocity_field(X, Y, params)

    # Grid spacing
    dx = (params['x_max'] - params['x_min']) / (params['grid_nx'] - 1)
    dy = (params['y_max'] - params['y_min']) / (params['grid_ny'] - 1)

    # ∂v/∂x using central differences
    dvdx = np.zeros_like(V)
    dvdx[:, 1:-1] = (V[:, 2:] - V[:, :-2]) / (2 * dx)
    dvdx[:, 0] = (V[:, 1] - V[:, 0]) / dx      # Forward at boundary
    dvdx[:, -1] = (V[:, -1] - V[:, -2]) / dx   # Backward at boundary

    # ∂u/∂y using central differences
    dudy = np.zeros_like(U)
    dudy[1:-1, :] = (U[2:, :] - U[:-2, :]) / (2 * dy)
    dudy[0, :] = (U[1, :] - U[0, :]) / dy
    dudy[-1, :] = (U[-1, :] - U[-2, :]) / dy

    omega = dvdx - dudy
    return omega


# ============================================================================
# DECAY MODEL
# ============================================================================

def compute_decay(x, params):
    """
    Model exponential decay of vortex street downstream.

    Physics: Viscous dissipation causes vorticity decay as exp(-α·x)
    where α ∝ ν/(U·D²) ∝ 1/Re

    Empirical: Decay length L_decay ~ Re·D (higher Re → longer street)
    """
    Re = params['Re']
    alpha = params['decay_alpha'] / Re  # Decay rate inversely proportional to Re
    x0 = params['decay_x0']

    # Exponential decay starting at x0
    decay_envelope = np.exp(-alpha * np.maximum(0, x - x0))

    # Decay length: distance for amplitude to drop to 1/e
    L_decay = 1.0 / alpha if alpha > 0 else np.inf

    return decay_envelope, L_decay


# ============================================================================
# TIME-VARYING DYNAMICS
# ============================================================================

def simulate_dynamic_vortex_street(params):
    """
    Main time-stepping loop for dynamic vortex tracking.

    Physics: Vortices are shed periodically at Strouhal frequency, then advected
    by the velocity field using RK4 integration. Cores grow via Lamb-Oseen diffusion.

    Returns:
        times: np.array of time values [T = D/U∞]
        vortices_history: List of vortex lists at each timestep
        u_probe: np.array of u velocity at probe
        v_probe: np.array of v velocity at probe
    """
    dt = params['dt']
    t_max = params['t_max']
    n_steps = int(t_max / dt)

    x_probe = params['probe_x']
    y_probe = params['probe_y']
    x_max = params['x_max']
    y_max = params['y_max']

    # Storage
    times = []
    vortices_history = []
    u_probe_history = []
    v_probe_history = []

    # Initialize
    vortices = []
    shedding_state = {'count': 0}
    t = 0.0

    print(f"\n{'='*70}")
    print(f"DYNAMIC VORTEX SIMULATION")
    print(f"{'='*70}")
    print(f"Steps: {n_steps}, dt={dt:.4f}, t_max={t_max}")
    print(f"Strouhal: St={params['St']:.3f}")
    print(f"Circulation coeff: C_gamma={params['C_gamma']:.3f}")
    print(f"Shedding period: T={params['D']/(params['St']*params['U_inf']):.3f}")
    print(f"Image vortices: {params.get('use_image_vortices', False)}")
    print(f"{'='*70}\n")

    # Time loop
    for step in range(n_steps):
        # 1. Shed vortex if needed
        shed_vortex(vortices, t, params, shedding_state)

        # 2. Advect vortices (RK4)
        if len(vortices) > 0:
            advect_vortices_rk4(vortices, t, dt, params)

        # 3. Remove vortices outside domain (for performance)
        vortices = [v for v in vortices
                    if v.x < x_max and abs(v.y) < y_max]

        # 4. Sample velocity at probe point
        u_p = params['U_inf']
        v_p = 0.0

        for vortex in vortices:
            # Direct contribution
            dx = x_probe - vortex.x
            dy = y_probe - vortex.y
            r_sq = dx**2 + dy**2 + 1e-10
            r = np.sqrt(r_sq)
            sigma = vortex.sigma(t, params)
            exp_factor = np.exp(-r_sq / sigma**2)
            vel_mag = (vortex.gamma / (2 * np.pi * r)) * (1 - exp_factor)
            u_p += vel_mag * (-dy / r)
            v_p += vel_mag * (dx / r)

            # Image contribution
            if params.get('use_image_vortices', False):
                image_pos = get_image_vortex_position(vortex, params)
                if image_pos is not None:
                    x_img, y_img = image_pos
                    dx_img = x_probe - x_img
                    dy_img = y_probe - y_img
                    r_sq_img = dx_img**2 + dy_img**2 + 1e-10
                    r_img = np.sqrt(r_sq_img)
                    vel_mag_img = (-vortex.gamma / (2 * np.pi * r_img)) * (1 - np.exp(-r_sq_img / sigma**2))
                    u_p += vel_mag_img * (-dy_img / r_img)
                    v_p += vel_mag_img * (dx_img / r_img)

        # 5. Store snapshot (deep copy vortices)
        times.append(t)
        vortices_history.append([Vortex(v.x, v.y, v.gamma, v.t_birth)
                                for v in vortices])
        u_probe_history.append(u_p)
        v_probe_history.append(v_p)

        # 6. Progress output
        if step % 500 == 0 or step == n_steps - 1:
            print(f"Step {step:5d}/{n_steps}: t={t:7.2f}, N_vortices={len(vortices):4d}, "
                  f"u_probe={u_p:+.3f}, v_probe={v_p:+.3f}")

        # 7. Advance time
        t += dt

    print(f"\n{'='*70}")
    print(f"Simulation complete: {len(vortices_history)} timesteps saved")
    print(f"Final vortex count: {len(vortices)}")
    print(f"{'='*70}\n")

    return (np.array(times), vortices_history,
            np.array(u_probe_history), np.array(v_probe_history))


# ============================================================================
# VISUALIZATION
# ============================================================================

def create_visualization(params):
    """
    Generate 4 separate visualizations of vortex street.

    Plots:
    [1] Vorticity field with streamlines, markers, cone, cylinder
    [2] Velocity magnitude contours
    [3] Downstream decay with Reynolds number sweep
    [4] Time-varying velocity at probe
    """
    # Create meshgrid
    x = np.linspace(params['x_min'], params['x_max'], params['grid_nx'])
    y = np.linspace(params['y_min'], params['y_max'], params['grid_ny'])
    X, Y = np.meshgrid(x, y)

    # Compute fields for base case
    U, V = compute_velocity_field(X, Y, params)
    omega = compute_vorticity_field(X, Y, params)
    vel_mag = np.sqrt(U**2 + V**2)
    vortex_pos = get_vortex_positions(params)

    # List to store all figures
    figures = []

    # ========================================================================
    # PLOT 1: Vorticity + Streamlines + Geometry
    # ========================================================================
    fig1 = plt.figure(figsize=(12, 8))
    ax1 = fig1.add_subplot(111)

    # Vorticity contours with diverging colormap
    vort_max = np.max(np.abs(omega))
    vort_levels = np.linspace(-vort_max, vort_max, 25)
    cf1 = ax1.contourf(X, Y, omega, levels=vort_levels, cmap='RdBu_r', extend='both')
    plt.colorbar(cf1, ax=ax1, label='Vorticity ω [1/s]')

    # Streamlines
    ax1.streamplot(X, Y, U, V, color='black', linewidth=0.5, density=1.5, arrowsize=0.8)

    # Vortex position markers
    ax1.scatter(vortex_pos['x_pos'], vortex_pos['y_pos'], c='red', marker='+',
                s=60, linewidths=2, label='+ Circulation', zorder=5)
    ax1.scatter(vortex_pos['x_neg'], vortex_pos['y_neg'], c='blue', marker='_',
                s=60, linewidths=2, label='- Circulation', zorder=5)

    # Cone boundaries
    x_cone = np.linspace(0, params['x_max'], 100)
    y_upper, y_lower = get_cone_boundaries(params, x_cone)
    ax1.plot(x_cone, y_upper, 'g--', linewidth=2, label='Cone boundary', zorder=4)
    ax1.plot(x_cone, y_lower, 'g--', linewidth=2, zorder=4)

    # Cylinder at origin
    cylinder = Circle((0, 0), params['D']/2, color='gray', fill=True, zorder=10)
    ax1.add_patch(cylinder)

    ax1.set_xlabel('x/D', fontsize=14)
    ax1.set_ylabel('y/D', fontsize=14)
    ax1.set_title('Vorticity Field with Streamlines', fontsize=16, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=12)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    plt.tight_layout()
    figures.append(('vorticity_field.png', fig1))

    # ========================================================================
    # PLOT 2: Velocity Magnitude (Focused on relevant region)
    # ========================================================================
    fig2 = plt.figure(figsize=(14, 6))
    ax2 = fig2.add_subplot(111)

    # Calculate decay length to determine plot limits
    _, L_decay = compute_decay(0, params)
    x_max_plot = min(params['x_max'], params['decay_x0'] + 2 * L_decay)

    # Create focused grid for detailed visualization
    x_focused = np.linspace(0, x_max_plot, 500)
    y_focused = np.linspace(params['y_min'], params['y_max'], 200)
    X_focused, Y_focused = np.meshgrid(x_focused, y_focused)

    # Compute velocity field on focused grid
    U_focused, V_focused = compute_velocity_field(X_focused, Y_focused, params)
    vel_mag_focused = np.sqrt(U_focused**2 + V_focused**2)

    # Add shaded region for 0 < x/D < 5 (model not suitable for near wake)
    ax2.axvspan(0, 5.0, alpha=0.12, color='red', zorder=1,
                label='Near wake (model less accurate, 0<x/D<5)')

    # Add vertical line at decay onset
    ax2.axvline(params['decay_x0'], color='orange', linestyle='--',
                linewidth=2, alpha=0.6, label=f'Decay onset (x/D={params["decay_x0"]})')

    # Velocity magnitude contours with more levels for detail
    cf2 = ax2.contourf(X_focused, Y_focused, vel_mag_focused, levels=50, cmap='viridis')
    cbar = plt.colorbar(cf2, ax=ax2, label='Velocity Magnitude |v| [m/s]')

    # Plot DIVERGING vortex street boundaries (connecting actual vortex positions)
    # This shows the REALISTIC wake spreading, not parallel rows!
    x_boundary = np.linspace(0, x_max_plot, 200)
    y_upper_boundary = np.array([compute_transverse_spacing(x, params) / 2 for x in x_boundary])
    y_lower_boundary = np.array([-compute_transverse_spacing(x, params) / 2 for x in x_boundary])
    ax2.plot(x_boundary, y_upper_boundary, 'lime', linewidth=2, linestyle='--',
             alpha=0.9, label='Diverging wake boundaries (viscous model)')
    ax2.plot(x_boundary, y_lower_boundary, 'lime', linewidth=2, linestyle='--', alpha=0.9)

    # Vortex markers (only for vortices in visible range) - DEFINE FIRST!
    visible_pos = vortex_pos['x_pos'][(vortex_pos['x_pos'] >= 0) & (vortex_pos['x_pos'] <= x_max_plot)]
    visible_pos_y = vortex_pos['y_pos'][(vortex_pos['x_pos'] >= 0) & (vortex_pos['x_pos'] <= x_max_plot)]
    visible_neg = vortex_pos['x_neg'][(vortex_pos['x_neg'] >= 0) & (vortex_pos['x_neg'] <= x_max_plot)]
    visible_neg_y = vortex_pos['y_neg'][(vortex_pos['x_neg'] >= 0) & (vortex_pos['x_neg'] <= x_max_plot)]

    # Optional: Draw vortex core size circles (shows growing cores downstream)
    # Plot every 3rd vortex to avoid clutter
    for i in range(0, len(visible_pos), 3):
        if i < len(visible_pos):
            x_v, y_v = visible_pos[i], visible_pos_y[i]
            sigma_v = compute_vortex_core_radius(x_v, params)
            core_circle = Circle((x_v, y_v), sigma_v, color='red', fill=False,
                               linestyle=':', linewidth=0.8, alpha=0.5)
            ax2.add_patch(core_circle)
    for i in range(0, len(visible_neg), 3):
        if i < len(visible_neg):
            x_v, y_v = visible_neg[i], visible_neg_y[i]
            sigma_v = compute_vortex_core_radius(x_v, params)
            core_circle = Circle((x_v, y_v), sigma_v, color='blue', fill=False,
                               linestyle=':', linewidth=0.8, alpha=0.5)
            ax2.add_patch(core_circle)

    ax2.scatter(visible_pos, visible_pos_y, c='red', marker='+',
                s=50, linewidths=2, zorder=5)
    ax2.scatter(visible_neg, visible_neg_y, c='blue', marker='_',
                s=50, linewidths=2, zorder=5)

    # Cylinder
    cylinder2 = Circle((0, 0), params['D']/2, color='gray', fill=True, zorder=10,
                       edgecolor='black', linewidth=2)
    ax2.add_patch(cylinder2)

    ax2.set_xlabel('x/D', fontsize=14, fontweight='bold')
    ax2.set_ylabel('y/D', fontsize=14, fontweight='bold')
    ax2.set_title('Velocity Magnitude - Viscous Gaussian Vortex Blob Model (Lamb-Oseen)',
                  fontsize=15, fontweight='bold')
    ax2.set_xlim([0, x_max_plot])
    ax2.set_ylim([params['y_min'], params['y_max']])
    ax2.legend(loc='upper right', fontsize=10, framealpha=0.9)
    ax2.grid(True, alpha=0.3, linewidth=0.5)
    plt.tight_layout()
    figures.append(('velocity_magnitude.png', fig2))

    # ========================================================================
    # PLOT 3: Downstream Decay with Reynolds Sweep
    # ========================================================================
    fig3 = plt.figure(figsize=(12, 8))
    ax3 = fig3.add_subplot(111)

    # Reynolds numbers to sweep
    Re_values = [60, 100, 140, 180]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    # Store base case Re
    Re_base = params['Re']

    for Re, color in zip(Re_values, colors):
        # Update Reynolds number
        params['Re'] = Re

        # Recompute vorticity for this Re
        omega_Re = compute_vorticity_field(X, Y, params)

        # Extract max vorticity along x-axis
        omega_max_vs_x = np.max(np.abs(omega_Re), axis=0)

        # Compute decay envelope
        decay_env, L_decay = compute_decay(x, params)

        # Plot with label
        linewidth = 3 if Re == Re_base else 2
        linestyle = '-' if Re == Re_base else '--'
        label = f'Re={Re:.0f} (L={L_decay:.1f}D)'
        if Re == Re_base:
            label += ' [base]'

        ax3.plot(x, omega_max_vs_x, color=color, linewidth=linewidth,
                linestyle=linestyle, label=label, alpha=0.8)

    # Restore base Re
    params['Re'] = Re_base

    ax3.set_xlabel('Downstream Distance x/D', fontsize=14)
    ax3.set_ylabel('Max |ω| [1/s]', fontsize=14)
    ax3.set_title('Downstream Decay - Reynolds Number Sweep', fontsize=16, fontweight='bold')
    ax3.legend(loc='upper right', fontsize=12)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim([0, params['x_max']])
    ax3.set_ylim(bottom=0)
    plt.tight_layout()
    figures.append(('downstream_decay.png', fig3))

    # ========================================================================
    # PLOT 4: Time-Varying Velocity (Dynamic Simulation)
    # ========================================================================
    fig4 = plt.figure(figsize=(12, 8))
    ax4 = fig4.add_subplot(111)

    print("\nRunning dynamic vortex simulation...")
    times, vortices_history, u_probe, v_probe = simulate_dynamic_vortex_street(params)
    vel_mag_probe = np.sqrt(u_probe**2 + v_probe**2)

    # Plot velocity components
    ax4.plot(times, u_probe, 'b-', label='u (streamwise)', linewidth=2)
    ax4.plot(times, v_probe, 'r-', label='v (transverse)', linewidth=2)
    ax4.plot(times, vel_mag_probe, 'g--', label='|v| (magnitude)',
             linewidth=1.5, alpha=0.7)

    # Detect shedding period from peaks
    from scipy import signal
    peaks, properties = signal.find_peaks(u_probe, prominence=0.1*np.std(u_probe))
    St_measured = np.nan
    T_measured = np.nan

    if len(peaks) > 1:
        T_measured = np.mean(np.diff(times[peaks]))
        f_measured = 1.0 / T_measured
        St_measured = f_measured * params['D'] / params['U_inf']

        # Mark first few periods
        for i in range(min(3, len(peaks))):
            ax4.axvline(times[peaks[i]], color='black', linestyle=':', alpha=0.3)

        label_text = f'Measured: T={T_measured:.3f}, St={St_measured:.3f}'
        ax4.axvline(times[peaks[0]], color='black', linestyle=':',
                    alpha=0.5, label=label_text)

    # Theoretical period
    T_theory = params['D'] / (params['St'] * params['U_inf'])
    St_theory = params['St']

    ax4.set_xlabel('Time t [D/U∞]', fontsize=14)
    ax4.set_ylabel('Velocity [m/s]', fontsize=14)
    ax4.set_title(f'Dynamic Simulation: Probe at ({params["probe_x"]}D, {params["probe_y"]}D)\n'
                  f'Input: St={St_theory:.3f}, T={T_theory:.3f} | '
                  f'Measured: St={St_measured:.3f}, T={T_measured:.3f}',
                  fontsize=13, fontweight='bold')
    ax4.legend(loc='upper right', fontsize=11)
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim([0, min(50, params['t_max'])])  # Show first 50 time units
    plt.tight_layout()
    figures.append(('time_varying_velocity_dynamic.png', fig4))

    # ========================================================================
    # Save Velocity Field Snapshots to Folder
    # ========================================================================
    import os

    # Create snapshots folder
    snapshot_folder = params['snapshots_folder']
    os.makedirs(snapshot_folder, exist_ok=True)

    # Select evenly-spaced snapshots
    n_snapshots = params['n_snapshots']
    snapshot_indices = np.linspace(0, len(times)-1, n_snapshots, dtype=int)

    print(f"\nSaving {n_snapshots} velocity field snapshots to '{snapshot_folder}/'...")

    # Create meshgrid for velocity field
    x_snap = np.linspace(params['x_min'], params['x_max'], params['grid_nx'])
    y_snap = np.linspace(params['y_min'], params['y_max'], params['grid_ny'])
    X_snap, Y_snap = np.meshgrid(x_snap, y_snap)

    for idx, snap_idx in enumerate(snapshot_indices):
        t_snap = times[snap_idx]
        vortices_snap = vortices_history[snap_idx]

        # Compute velocity field at this snapshot
        U_snap = np.full_like(X_snap, params['U_inf'])
        V_snap = np.zeros_like(X_snap)

        for vortex in vortices_snap:
            U_v, V_v = gaussian_velocity_contribution_dynamic(X_snap, Y_snap, vortex, t_snap, params)
            U_snap += U_v
            V_snap += V_v

            # Add image contribution if enabled
            if params.get('use_image_vortices', False):
                image_pos = get_image_vortex_position(vortex, params)
                if image_pos is not None:
                    x_img, y_img = image_pos
                    sigma_img = vortex.sigma(t_snap, params)
                    U_img, V_img = gaussian_velocity_contribution(
                        X_snap, Y_snap, x_img, y_img, sigma_img, -vortex.gamma
                    )
                    U_snap += U_img
                    V_snap += V_img

        # Create individual snapshot figure
        # Wider figure to accommodate equal aspect ratio and match colorbar height
        fig_snap = plt.figure(figsize=(18, 6))
        ax_snap = fig_snap.add_subplot(111)

        # Velocity magnitude contour
        vel_mag = np.sqrt(U_snap**2 + V_snap**2)
        cf = ax_snap.contourf(X_snap, Y_snap, vel_mag, levels=50, cmap='viridis')
        plt.colorbar(cf, ax=ax_snap, label='Velocity Magnitude |v| [m/s]')

        # Velocity vector field (quiver with reduced density and smaller arrows)
        skip = 8  # Plot every 8th vector (reduced density)
        ax_snap.quiver(X_snap[::skip, ::skip], Y_snap[::skip, ::skip],
                      U_snap[::skip, ::skip], V_snap[::skip, ::skip],
                      color='white', alpha=0.6, scale=35, width=0.002)

        # Cylinder
        cylinder = Circle((0, 0), params['D']/2, color='gray',
                         fill=True, zorder=10, edgecolor='black', linewidth=2)
        ax_snap.add_patch(cylinder)

        # Formatting
        ax_snap.set_xlim([params['x_min'], params['x_max']])
        ax_snap.set_ylim([params['y_min'], params['y_max']])
        ax_snap.set_xlabel('x/D', fontsize=12, fontweight='bold')
        ax_snap.set_ylabel('y/D', fontsize=12, fontweight='bold')
        ax_snap.set_title(f'Velocity Field at t = {t_snap:.2f} [D/U∞] '
                         f'(snapshot {idx+1}/{n_snapshots}, N_vortices={len(vortices_snap)})',
                         fontsize=13, fontweight='bold')
        ax_snap.set_aspect('equal')
        ax_snap.grid(True, alpha=0.3)

        # Save to folder
        filename = os.path.join(snapshot_folder, f'velocity_field_t{t_snap:06.2f}.png')
        fig_snap.savefig(filename, dpi=200, bbox_inches='tight')
        plt.close(fig_snap)

        print(f"  Saved: {filename}")

    print(f"All {n_snapshots} snapshots saved to '{snapshot_folder}/'")

    return figures


# ============================================================================
# METRICS AND VALIDATION
# ============================================================================

def compute_metrics(params):
    """
    Calculate and display key physics metrics.

    Validation:
    - V₀ from orbital velocity formula
    - Strouhal number St = f·D/U (expect ~0.2 for Re=100)
    - Cone angle from geometry
    - Decay length from exponential model
    """
    h = params['h']
    l = params['l']
    Gamma = params['Gamma']
    D = params['D']
    U_inf = params['U_inf']
    Re = params['Re']

    # Orbital velocity (tanh formula for a = l/2)
    V0 = (Gamma / (2 * l)) * np.tanh(np.pi * h / l)

    # Shedding frequency and Strouhal number
    f_shed = V0 / l
    St = f_shed * D / U_inf

    # Cone angle
    theta_deg = compute_cone_angle(params)['theta_full_deg']

    # Decay length
    _, L_decay = compute_decay(0, params)

    print("=" * 70)
    print("VISCOUS GAUSSIAN VORTEX BLOB MODEL - PHYSICS METRICS")
    print("(Lamb-Oseen vortex with wake spreading)")
    print("=" * 70)
    print(f"Orbital Velocity V0:      {V0:.4f} m/s")
    print(f"Shedding Frequency f:     {f_shed:.4f} Hz")
    print(f"Strouhal Number St:       {St:.4f} (expected ~0.2 for Re={Re:.0f})")
    print(f"Cone Angle (full):        {theta_deg:.2f} degrees")
    print(f"Decay Length:             {L_decay:.2f} D")
    print(f"h/l Spacing Ratio:        {h/l:.4f} (stability: 0.281)")
    print(f"Reynolds Number:          {Re:.1f}")
    print()
    print("Viscous Model Parameters:")
    print(f"  Initial core radius:    {params['sigma_0']:.3f} D")
    print(f"  Spreading coefficient:  {params['spreading_coef']:.3f}")
    print(f"  Core at x=5D:           {compute_vortex_core_radius(5.0, params):.3f} D")
    print(f"  Wake width at x=10D:    {compute_transverse_spacing(10.0, params):.3f} D")
    print("=" * 70)

    return {
        'V0': V0,
        'frequency': f_shed,
        'Strouhal': St,
        'cone_angle': theta_deg,
        'decay_length': L_decay
    }


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """
    Main execution: compute metrics, generate visualization, save output.
    """
    # Compute and display physics metrics
    metrics = compute_metrics(params)

    # Generate 4 separate visualizations
    print("\nGenerating visualizations...")
    figures = create_visualization(params)

    # Save all high-resolution figures
    print("\nSaving plots:")
    for filename, fig in figures:
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Saved: {filename} (300 dpi)")

    # Display interactive plots
    plt.show()


if __name__ == '__main__':
    main()
