"""Viscous vortex method for Von Kármán vortex streets shed from one or more cylinders."""

import os
from datetime import datetime

import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.interpolate import interp1d
from numba import njit, prange

# CONSTANTS
INITIAL_CORE_SIZE_FACTOR = 0.1
CIRCULATION_FACTOR = 1.5
CYLINDER_TOLERANCE = 0.9
SIGMOID_STEEPNESS = 10.0
FORMATION_NUMBER = 4.0
WAKE_TRANSVERSE_BASE = 0.5
WAKE_TRANSVERSE_GROWTH = 0.1
WAKE_STREAMWISE_DECAY = 20.0
THETA_SEP_STD_DEG = 5.0
HIGH_RE_FORMATION_REDUCTION = 0.8

# PHYSICS FUNCTIONS
def initial_core_size(D, Re):
    """Compute initial vortex core size at shedding."""
    return INITIAL_CORE_SIZE_FACTOR * D


def compute_strouhal_number(Re):
    """Strouhal number from Reynolds (Roshko 1961 correlations)."""
    if Re < 47:
        return 0.0
    elif Re < 300:
        return 0.212 * (1 - 21.2 / Re)
    elif Re < 2e5:
        return 0.20
    elif Re < 3.5e6:
        log_Re = np.log10(Re)
        return 0.20 + 0.06 * (log_Re - 5.3)
    else:
        return 0.30


def compute_circulation_magnitude(U_inf, D):
    """Target circulation Gamma."""
    return CIRCULATION_FACTOR * U_inf * D


def compute_shedding_period(U_inf, D, St):
    """Shedding period T = D/(St*U) from Strouhal relation f = St*U/D."""
    if St <= 0:
        return np.inf
    return D / (St * U_inf)


def compute_eddy_viscosity_field(x, y, D, nu_molecular, U_inf, Re):
    """Effective viscosity (molecular + turbulent eddy contribution in wake)."""
    if Re < 1000:
        C_eddy = 0.01
    elif Re < 1e5:
        C_eddy = 0.02 + 0.01 * np.log10(Re / 1000)
    else:
        C_eddy = 0.10

    x_norm = x / D
    y_norm = y / D

    sigma_wake = WAKE_TRANSVERSE_BASE + WAKE_TRANSVERSE_GROWTH * x_norm
    f_transverse = np.exp(-(y_norm**2) / (2 * sigma_wake**2))
    f_streamwise = np.exp(-x_norm / WAKE_STREAMWISE_DECAY)
    f_wake = f_transverse * f_streamwise

    nu_turbulent = C_eddy * U_inf * D * f_wake
    return nu_molecular + nu_turbulent


def compute_formation_time(D, U_inf, Re):
    """Vortex formation time scale (formation number F* ~ 4.0)."""
    F_star = FORMATION_NUMBER
    if Re > 1e5:
        F_star *= HIGH_RE_FORMATION_REDUCTION
    return F_star * (D / U_inf)


def circulation_growth_sigmoid(t_age, T_form):
    """Smooth circulation growth factor [0,1] during vortex formation."""
    if t_age < 0:
        return 0.0

    tau = t_age / T_form
    k = SIGMOID_STEEPNESS

    sigmoid_raw = 1.0 / (1.0 + np.exp(-k * (tau - 0.5)))
    g_min = 1.0 / (1.0 + np.exp(k * 0.5))
    g_max = 1.0 / (1.0 + np.exp(-k * 0.5))
    g_normalized = (sigmoid_raw - g_min) / (g_max - g_min)

    return np.clip(g_normalized, 0.0, 1.0)



# NUMBA KERNEL
@njit(parallel=True, fastmath=True)
def compute_velocities_numba(targets_x, targets_y, sources_x, sources_y,
                             gammas, sigmas, U_inf_x, U_inf_y):
    """Numba kernel: velocity induced by source vortices at target points."""
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

            if R_sq > 1e-12:
                R = np.sqrt(R_sq)
                sigma_sq = sigmas[i] * sigmas[i]
                EXP = np.exp(-R_sq / sigma_sq)
                V_theta = (gammas[i] / (2.0 * np.pi * R)) * (1.0 - EXP)

                u_total += -V_theta * DY / R
                v_total += V_theta * DX / R

        U[j] = u_total
        V[j] = v_total

    return U, V



# HELPER FUNCTIONS
def apply_transforms(x, y, rot_deg, flow_deg):
    """Apply rotation and flow angle transforms to coordinates."""
    total_angle = np.radians(rot_deg + flow_deg)
    c, s = np.cos(total_angle), np.sin(total_angle)
    return c*x - s*y, s*x + c*y


def _freestream_components(U_inf, flow_angle):
    """Decompose freestream into (Ux, Uy) from angle in degrees."""
    rad = np.radians(flow_angle)
    return U_inf * np.cos(rad), U_inf * np.sin(rad)


def compute_velocities_at_points(targets_x, targets_y, sources_x, sources_y,
                                 gammas, sigmas, U_inf, flow_angle):
    """Calculate velocity at arbitrary target points due to source vortices."""
    if len(targets_x) == 0:
        return np.array([]), np.array([])
    U_inf_x, U_inf_y = _freestream_components(U_inf, flow_angle)
    return compute_velocities_numba(
        targets_x, targets_y, sources_x, sources_y,
        gammas, sigmas, U_inf_x, U_inf_y
    )


def compute_velocity_field(X, Y, vortices, U_inf, flow_angle):
    """Compute velocity field over a meshgrid.

    ``vortices`` is a list of dicts with keys 'x', 'y', 'gamma', 'sigma'
    (the snapshot format stored in results_df).
    """
    vortex_x = np.array([v['x'] for v in vortices])
    vortex_y = np.array([v['y'] for v in vortices])
    vortex_gamma = np.array([v['gamma'] for v in vortices])
    vortex_sigma = np.array([v['sigma'] for v in vortices])

    U_inf_x, U_inf_y = _freestream_components(U_inf, flow_angle)
    X_flat = X.ravel()
    Y_flat = Y.ravel()

    U_flat, V_flat = compute_velocities_numba(
        X_flat, Y_flat, vortex_x, vortex_y,
        vortex_gamma, vortex_sigma, U_inf_x, U_inf_y
    )
    return U_flat.reshape(X.shape), V_flat.reshape(X.shape)


def shed_vortex_with_turbulence(cyl, t, upper, Re, U_inf, D_ref, Gamma_mag,
                                rotation_angle, flow_angle, theta_sep,
                                enable_stochastic=True):
    """Create new vortex at separation point with optional turbulent fluctuations."""
    sign = 1 if upper else -1
    gamma_base = -Gamma_mag if upper else Gamma_mag

    if enable_stochastic and Re > 1000:
        fluctuation_level = min(0.10 * np.log10(Re / 1000), 0.20)
        gamma = gamma_base * (1 + np.random.normal(0, fluctuation_level))
        theta_sep = theta_sep + np.random.normal(0, np.radians(THETA_SEP_STD_DEG))
    else:
        gamma = gamma_base

    x_cyl_global, y_cyl_global = apply_transforms(
        cyl['x'], cyl['y'], rotation_angle, flow_angle
    )

    rear_angle = np.radians(flow_angle) + np.pi
    theta = rear_angle - sign * theta_sep
    a = cyl['D'] / 2

    return {
        'x': x_cyl_global + a * np.cos(theta),
        'y': y_cyl_global + a * np.sin(theta),
        'gamma_target': gamma,
        'gamma': 0.0,
        'sigma': initial_core_size(D_ref, Re),
        'sigma_0': initial_core_size(D_ref, Re),
        'birth_t': t,
        'U_inf_birth': U_inf,
        'T_form': compute_formation_time(D_ref, U_inf, Re),
    }


def remove_internal_vortices(cylinders, vortex_state, rotation_angle, flow_angle):
    """Remove vortices that have drifted inside any cylinder."""
    keep = np.ones(len(vortex_state), dtype=bool)

    for cyl in cylinders:
        x_cyl, y_cyl = apply_transforms(
            cyl['x'], cyl['y'], rotation_angle, flow_angle
        )
        dx = vortex_state.x - x_cyl
        dy = vortex_state.y - y_cyl
        dist = np.sqrt(dx**2 + dy**2)
        keep &= dist >= (cyl['D'] / 2) * CYLINDER_TOLERANCE

    vortex_state.filter(keep)



# VORTEX STATE (parallel NumPy arrays)
class VortexState:
    """Manages vortex data as parallel numpy arrays."""

    _FIELDS = ('x', 'y', 'gamma', 'gamma_target', 'sigma', 'sigma_0',
               'birth_t', 'U_inf_birth', 'T_form')

    def __init__(self):
        for f in self._FIELDS:
            setattr(self, f, np.empty(0, dtype=np.float64))

    def __len__(self):
        return len(self.x)

    def append(self, vortex_dict):
        """Append a single vortex (dict with matching keys)."""
        for f in self._FIELDS:
            current = getattr(self, f)
            setattr(self, f, np.append(current, vortex_dict[f]))

    def filter(self, mask):
        """Keep only vortices where *mask* is True."""
        for f in self._FIELDS:
            setattr(self, f, getattr(self, f)[mask])

    def snapshot(self, t):
        """Return list-of-dicts snapshot (for storage in results_df)."""
        return [
            {'x': self.x[i], 'y': self.y[i], 'gamma': self.gamma[i],
             'sigma': self.sigma[i], 'birth_t': self.birth_t[i],
             'age': t - self.birth_t[i]}
            for i in range(len(self))
        ]



# VORTEX SIMULATION CLASS
class VortexAmp:
    """Von Kármán vortex street simulation using viscous vortex method."""

    def __init__(self, cylinders, nu=4.88e-6, flow_angle_metocean=0,
                 rotation_angle=0, theta_sep_deg=80.0, x_removal=250.0,
                 dt=0.01, measurement_points=None, save_interval=0.0,
                 turbulence_thresholds=None, sigma_max_factor=0.5,
                 output_base_dir='output'):

        self.cylinders = cylinders
        self.output_base_dir = output_base_dir
        self.output_dir = None
        self.nu = nu
        self.flow_angle_metocean = flow_angle_metocean
        self.rotation_angle = rotation_angle
        self.theta_sep_deg = theta_sep_deg
        self.theta_sep = np.radians(theta_sep_deg)
        self.x_removal = x_removal
        self.dt = dt
        self.measurement_points = measurement_points or []
        self.save_interval = save_interval
        self.sigma_max_factor = sigma_max_factor

        self.D_ref = max(cyl['D'] for cyl in cylinders)
        self.flow_angle = 270.0 - flow_angle_metocean

        if turbulence_thresholds is None:
            turbulence_thresholds = {
                'eddy_viscosity': 10_000_000,
                'stochastic_shedding': 10_000_000,
                'core_saturation': 10_000_000,
            }
        self.turbulence_thresholds = turbulence_thresholds

        self.results_df = None

    # Velocity / flow-angle profile setup
    @staticmethod
    def _setup_velocity_profile(velocity_mode, **kwargs):
        """Return a callable ``get_velocity(t) -> float``."""
        if velocity_mode == 'constant':
            U_inf = kwargs.get('U_inf', 1.0)
            return lambda t: U_inf

        if velocity_mode == 'file':
            velocity_file = kwargs.get('velocity_file')
            if velocity_file is None:
                raise ValueError("velocity_file must be provided for mode='file'")
            df = pd.read_excel(velocity_file)
            interp = interp1d(df.iloc[:, 0].values, df.iloc[:, 1].values,
                              kind='linear', fill_value='extrapolate')
            return lambda t: float(interp(t))

        if velocity_mode == 'function':
            fn = kwargs.get('velocity_function')
            if fn is None:
                raise ValueError("velocity_function must be provided for mode='function'")
            return fn

        raise ValueError(f"Invalid velocity_mode: {velocity_mode}")

    @staticmethod
    def _setup_flow_angle_profile(flow_angle_mode, flow_angle_metocean, **kwargs):
        """Return a callable ``get_flow_angle(t) -> float`` (internal angle)."""
        if flow_angle_mode == 'constant':
            angle = kwargs.get('flow_angle_metocean', flow_angle_metocean)
            internal = 270.0 - angle
            return lambda t: internal

        if flow_angle_mode == 'file':
            flow_angle_file = kwargs.get('flow_angle_file')
            if flow_angle_file is None:
                raise ValueError("flow_angle_file must be provided for mode='file'")
            df = pd.read_excel(flow_angle_file)
            interp = interp1d(df.iloc[:, 0].values, df.iloc[:, 1].values,
                              kind='linear', fill_value='extrapolate')
            return lambda t: 270.0 - float(interp(t))

        if flow_angle_mode == 'function':
            fn = kwargs.get('flow_angle_function')
            if fn is None:
                raise ValueError("flow_angle_function must be provided for mode='function'")
            return lambda t: 270.0 - fn(t)

        raise ValueError(f"Invalid flow_angle_mode: {flow_angle_mode}")

    # Simulation sub-steps
    def _shed_vortices(self, t, vortices, shed_counters, next_shed_times,
                       U_inf_now, Re_now, St_now, flow_angle_now,
                       enable_stochastic):
        """Shed new vortices from each cylinder when due."""
        Gamma_mag = compute_circulation_magnitude(U_inf_now, self.D_ref)

        for c_idx, cyl in enumerate(self.cylinders):
            if t < next_shed_times[c_idx]:
                continue

            upper = (shed_counters[c_idx] % 2 == 0)
            new_vortex = shed_vortex_with_turbulence(
                cyl, t, upper, Re_now, U_inf_now, self.D_ref, Gamma_mag,
                self.rotation_angle, flow_angle_now, self.theta_sep,
                enable_stochastic=enable_stochastic,
            )
            vortices.append(new_vortex)
            shed_counters[c_idx] += 1

            shed_period = compute_shedding_period(U_inf_now, self.D_ref, St_now)
            next_shed_times[c_idx] += shed_period / 2  # two vortices per period

    def _diffuse_cores(self, t, vortices, U_inf_now, Re_now,
                       enable_eddy, enable_saturation):
        """Update core sizes (diffusion) and circulation (growth)."""
        N = len(vortices)
        if N == 0:
            return

        ages = t - vortices.birth_t

        if enable_eddy:
            for i in range(N):
                nu_eff = compute_eddy_viscosity_field(
                    vortices.x[i], vortices.y[i],
                    self.D_ref, self.nu, U_inf_now, Re_now,
                )
                sigma_new = np.sqrt(vortices.sigma_0[i]**2 + 4 * nu_eff * ages[i])
                if enable_saturation:
                    sigma_new = min(sigma_new, self.sigma_max_factor * self.D_ref)
                vortices.sigma[i] = sigma_new
        else:
            vortices.sigma = np.sqrt(vortices.sigma_0**2 + 4 * self.nu * ages)

        for i in range(N):
            growth = circulation_growth_sigmoid(ages[i], vortices.T_form[i])
            vortices.gamma[i] = vortices.gamma_target[i] * growth

    def _advect_rk4(self, t, vortices, get_velocity, get_flow_angle):
        """Advance vortex positions using RK4 with proper intermediate sources."""
        N = len(vortices)
        if N == 0:
            return

        dt = self.dt
        sx, sy = vortices.x.copy(), vortices.y.copy()
        gammas, sigmas = vortices.gamma, vortices.sigma

        # Stage 1
        U_inf_0 = get_velocity(t)
        fa_0 = get_flow_angle(t)
        U1, V1 = compute_velocities_at_points(sx, sy, sx, sy, gammas, sigmas, U_inf_0, fa_0)
        k1x, k1y = dt * U1, dt * V1

        # Stage 2
        U_inf_mid = get_velocity(t + 0.5 * dt)
        fa_mid = get_flow_angle(t + 0.5 * dt)
        mid_x, mid_y = sx + 0.5 * k1x, sy + 0.5 * k1y
        U2, V2 = compute_velocities_at_points(mid_x, mid_y, mid_x, mid_y, gammas, sigmas, U_inf_mid, fa_mid)
        k2x, k2y = dt * U2, dt * V2

        # Stage 3
        mid_x, mid_y = sx + 0.5 * k2x, sy + 0.5 * k2y
        U3, V3 = compute_velocities_at_points(mid_x, mid_y, mid_x, mid_y, gammas, sigmas, U_inf_mid, fa_mid)
        k3x, k3y = dt * U3, dt * V3

        # Stage 4
        U_inf_end = get_velocity(t + dt)
        fa_end = get_flow_angle(t + dt)
        end_x, end_y = sx + k3x, sy + k3y
        U4, V4 = compute_velocities_at_points(end_x, end_y, end_x, end_y, gammas, sigmas, U_inf_end, fa_end)
        k4x, k4y = dt * U4, dt * V4

        vortices.x = sx + (k1x + 2*k2x + 2*k3x + k4x) / 6
        vortices.y = sy + (k1y + 2*k2y + 2*k3y + k4y) / 6

    def _remove_downstream(self, vortices):
        """Remove vortices past the downstream removal boundary."""
        vortices.filter(vortices.x <= self.x_removal)

    def _record_probes(self, vortices, U_inf_now, flow_angle_now, history):
        """Record velocities at measurement probe locations."""
        n_probes = len(self.measurement_points)
        if n_probes == 0:
            return

        if len(vortices) > 0:
            probe_x = np.array([p[0] for p in self.measurement_points])
            probe_y = np.array([p[1] for p in self.measurement_points])
            ux, uy = compute_velocities_at_points(
                probe_x, probe_y, vortices.x, vortices.y,
                vortices.gamma, vortices.sigma, U_inf_now, flow_angle_now,
            )
        else:
            Ux, Uy = _freestream_components(U_inf_now, flow_angle_now)
            ux = np.full(n_probes, Ux)
            uy = np.full(n_probes, Uy)

        for i in range(n_probes):
            history[f'probe_{i}_ux'].append(ux[i])
            history[f'probe_{i}_uy'].append(uy[i])
            history[f'probe_{i}_vmag'].append(np.sqrt(ux[i]**2 + uy[i]**2))

    # Main entry point
    def run(self, velocity_mode, total_time, flow_angle_mode='constant',
            output_dir=None, **kwargs):
        """Run simulation to completion."""
        # Output directory
        if output_dir is not None:
            self.output_dir = output_dir
        else:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            self.output_dir = os.path.join(self.output_base_dir, f'sim_{timestamp}')
        os.makedirs(self.output_dir, exist_ok=True)
        print(f"Output directory: {self.output_dir}")

        # Profile callables
        get_velocity = self._setup_velocity_profile(velocity_mode, **kwargs)
        get_flow_angle = self._setup_flow_angle_profile(
            flow_angle_mode, self.flow_angle_metocean, **kwargs
        )

        # Initialise state
        vortices = VortexState()
        t = 0.0
        shed_counters = [0] * len(self.cylinders)

        U_inf_0 = get_velocity(t)
        Re_0 = U_inf_0 * self.D_ref / self.nu
        St_0 = compute_strouhal_number(Re_0)
        shed_period_0 = compute_shedding_period(U_inf_0, self.D_ref, St_0)
        next_shed_times = [shed_period_0] * len(self.cylinders)

        self._print_header(U_inf_0, Re_0, St_0, shed_period_0)

        # History storage
        history = {
            'time': [], 'U_inf': [], 'Re': [], 'St': [],
            'n_vortices': [], 'vortex_field': [],
        }
        for i in range(len(self.measurement_points)):
            history[f'probe_{i}_ux'] = []
            history[f'probe_{i}_uy'] = []
            history[f'probe_{i}_vmag'] = []

        next_save_time = self.save_interval if self.save_interval > 0 else np.inf
        num_steps = int(total_time / self.dt)

        # Main loop
        for step in tqdm(range(num_steps), desc="Simulating", mininterval=0.5, unit="step"):
            U_inf_now = get_velocity(t)
            flow_angle_now = get_flow_angle(t)
            Re_now = U_inf_now * self.D_ref / self.nu
            St_now = compute_strouhal_number(Re_now)

            enable_eddy = Re_now > self.turbulence_thresholds['eddy_viscosity']
            enable_stochastic = Re_now > self.turbulence_thresholds['stochastic_shedding']
            enable_saturation = Re_now > self.turbulence_thresholds['core_saturation']

            self._shed_vortices(
                t, vortices, shed_counters, next_shed_times,
                U_inf_now, Re_now, St_now, flow_angle_now, enable_stochastic,
            )

            self._diffuse_cores(t, vortices, U_inf_now, Re_now,
                                enable_eddy, enable_saturation)

            self._advect_rk4(t, vortices, get_velocity, get_flow_angle)
            self._remove_downstream(vortices)

            # Record history
            history['time'].append(t)
            history['U_inf'].append(U_inf_now)
            history['Re'].append(Re_now)
            history['St'].append(St_now)
            history['n_vortices'].append(len(vortices))

            self._record_probes(vortices, U_inf_now, flow_angle_now, history)

            if t >= next_save_time or (self.save_interval == 0 and step == num_steps - 1):
                history['vortex_field'].append(vortices.snapshot(t))
                next_save_time += self.save_interval
            else:
                history['vortex_field'].append(None)

            t += self.dt

        print(f"Simulation complete. Total vortices in domain: {len(vortices)}")

        self.results_df = pd.DataFrame(history)
        return self.results_df

    @staticmethod
    def _print_header(U_inf_0, Re_0, St_0, shed_period_0):
        print("Starting vortex method simulation...")
        print(f"  Initial velocity U_inf(0) = {U_inf_0:.3f} m/s")
        print(f"  Initial Reynolds Re(0) = {Re_0:.2e}")
        print(f"  Strouhal number St(0) = {St_0:.3f}")
        print(f"  Initial shedding period = {shed_period_0:.3f} s\n")

    # Results saving/loading
    def save_results(self, filename):
        """Save results DataFrame to pickle file in the output directory."""
        if self.results_df is None:
            raise ValueError("No results to save. Run simulation first.")
        if self.output_dir is None:
            raise ValueError("No output directory set. Run simulation first.")
        filepath = os.path.join(self.output_dir, filename)
        self.results_df.to_pickle(filepath)
        print(f"Results saved to {filepath}")

    @staticmethod
    def load_results(filename):
        """Load results DataFrame from pickle file."""
        return pd.read_pickle(filename)
