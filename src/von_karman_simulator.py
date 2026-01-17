"""Viscous vortex method for Von Kármán vortex streets shed from one or more cylinders."""

import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.interpolate import interp1d
from numba import njit, prange


# PHYSICS FUNCTIONS

def initial_core_size(D, Re):
    """Compute initial vortex core size at shedding"""
    return 0.1 * D


def compute_strouhal_number(Re):
    """Strouhal number from Reynolds (Roshko 1961 correlations)"""
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
    if Re < 1000:
        C_eddy = 0.01
    elif Re < 1e5:
        C_eddy = 0.02 + 0.01 * np.log10(Re / 1000)
    else:
        C_eddy = 0.10

    x_norm = x / D
    y_norm = y / D

    sigma_wake = 0.5 + 0.1 * x_norm
    f_transverse = np.exp(-(y_norm**2) / (2 * sigma_wake**2))
    f_streamwise = np.exp(-x_norm / 20.0)
    f_wake = f_transverse * f_streamwise

    nu_turbulent = C_eddy * U_inf * D * f_wake
    return nu_molecular + nu_turbulent


def compute_formation_time(D, U_inf, Re):
    """Vortex formation time scale (formation number F* ~ 4.0)"""
    F_star = 4.0
    if Re > 1e5:
        F_star *= 0.8
    return F_star * (D / U_inf)


def circulation_growth_sigmoid(t_age, T_form):
    """Smooth circulation growth factor [0,1] during vortex formation"""
    if t_age < 0:
        return 0.0

    tau = t_age / T_form
    k = 10.0

    sigmoid_raw = 1.0 / (1.0 + np.exp(-k * (tau - 0.5)))

    g_min = 1.0 / (1.0 + np.exp(k * 0.5))
    g_max = 1.0 / (1.0 + np.exp(-k * 0.5))
    g_normalized = (sigmoid_raw - g_min) / (g_max - g_min)

    return np.clip(g_normalized, 0.0, 1.0)


# NUMBA KERNELS

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


# HELPER FUNCTIONS

def apply_transforms(x, y, rot_deg, flow_deg):
    """Apply rotation and flow angle transforms to coordinates"""
    total_angle = np.radians(rot_deg + flow_deg)
    c, s = np.cos(total_angle), np.sin(total_angle)
    return c*x - s*y, s*x + c*y


def shed_vortex_with_turbulence(cyl, t, upper, Re, U_inf, D_ref, Gamma_mag,
                                  rotation_angle, flow_angle, theta_sep, enable_stochastic=True):
    """Create new vortex at separation point with turbulent fluctuations"""
    sign = 1 if upper else -1
    gamma_base = -Gamma_mag if upper else Gamma_mag

    if enable_stochastic and Re > 1000:
        fluctuation_level = 0.10 * np.log10(Re / 1000)
        fluctuation_level = min(fluctuation_level, 0.20)
        gamma = gamma_base * (1 + np.random.normal(0, fluctuation_level))
    else:
        gamma = gamma_base

    if enable_stochastic and Re > 1000:
        theta_sep_perturbed = theta_sep + np.random.normal(0, np.radians(5))
    else:
        theta_sep_perturbed = theta_sep

    x_cyl_global, y_cyl_global = apply_transforms(cyl['x'], cyl['y'], rotation_angle, flow_angle)

    flow_direction_rad = np.radians(flow_angle)
    rear_angle = flow_direction_rad + np.pi
    theta = rear_angle - sign * theta_sep_perturbed

    a = cyl['D'] / 2
    x_shed = x_cyl_global + a * np.cos(theta)
    y_shed = y_cyl_global + a * np.sin(theta)

    sigma_0 = initial_core_size(D_ref, Re)

    return {
        'x': x_shed,
        'y': y_shed,
        'gamma_target': gamma,
        'gamma': 0.0,
        'sigma': sigma_0,
        'sigma_0': sigma_0,
        'birth_t': t,
        'U_inf_birth': U_inf,
        'T_form': compute_formation_time(D_ref, U_inf, Re)
    }


def remove_internal_vortices(cylinders, vortices, rotation_angle, flow_angle):
    """Filter vortices inside cylinders"""
    filtered_vortices = []
    tolerance = 0.9

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


def compute_velocities_at_points(targets_x, targets_y, sources_x, sources_y, gammas, sigmas, U_inf, flow_angle):
    """Calculate velocity at arbitrary points"""
    M = len(targets_x)
    if M == 0:
        return np.array([]), np.array([])

    flow_angle_rad = np.radians(flow_angle)
    U_inf_x = U_inf * np.cos(flow_angle_rad)
    U_inf_y = U_inf * np.sin(flow_angle_rad)

    return compute_velocities_numba(targets_x, targets_y, sources_x, sources_y, gammas, sigmas, U_inf_x, U_inf_y)


def compute_velocity_field(X, Y, vortices, U_inf, flow_angle):
    """Compute velocity field over meshgrid"""
    vortex_x = np.array([v['x'] for v in vortices])
    vortex_y = np.array([v['y'] for v in vortices])
    vortex_gamma = np.array([v['gamma'] for v in vortices])
    vortex_sigma = np.array([v['sigma'] for v in vortices])

    flow_angle_rad = np.radians(flow_angle)
    U_inf_x = U_inf * np.cos(flow_angle_rad)
    U_inf_y = U_inf * np.sin(flow_angle_rad)

    X_flat = X.ravel()
    Y_flat = Y.ravel()

    U_flat, V_flat = compute_velocity_field_numba(
        X_flat, Y_flat, vortex_x, vortex_y, vortex_gamma, vortex_sigma, U_inf_x, U_inf_y
    )

    U = U_flat.reshape(X.shape)
    V = V_flat.reshape(X.shape)

    return U, V


# VORTEX SIMULATION CLASS

class VortexAmp:
    """Von Kármán vortex street simulation using viscous vortex method"""

    def __init__(self, cylinders, nu=4.88e-6, flow_angle_metocean=0, rotation_angle=0,
                 theta_sep_deg=80.0, x_removal=250.0, dt=0.01,
                 measurement_points=None, save_interval=0.0,
                 turbulence_thresholds=None, sigma_max_factor=0.5):
        """
        Initialize vortex simulation.

        Parameters
        ----------
        cylinders : list[dict]
            List of cylinder definitions: [{'x': x, 'y': y, 'D': D}, ...]
        nu : float
            Kinematic viscosity [m²/s]
        flow_angle_metocean : float
            Flow direction [degrees, metocean convention]
        rotation_angle : float
            Cylinder cluster rotation [degrees]
        theta_sep_deg : float
            Separation angle [degrees] from rear stagnation point
        x_removal : float
            Remove vortices beyond this x-coordinate [m]
        dt : float
            Time step [s]
        measurement_points : list[tuple]
            List of (x, y) coordinates for velocity measurements
        save_interval : float
            Save vortex field every N seconds (0 = only save final state)
        turbulence_thresholds : dict
            Dict with 'eddy_viscosity', 'stochastic_shedding', 'core_saturation' Re thresholds
        sigma_max_factor : float
            Maximum core size = sigma_max_factor * D_ref
        """
        self.cylinders = cylinders
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
                'eddy_viscosity': 10000000,
                'stochastic_shedding': 10000000,
                'core_saturation': 10000000
            }
        self.turbulence_thresholds = turbulence_thresholds

        self.results_df = None

    def _setup_velocity_profile(self, velocity_mode, **velocity_kwargs):
        """Setup velocity profile based on mode"""
        if velocity_mode == 'constant':
            U_inf = velocity_kwargs.get('U_inf', 1.0)
            return lambda t: U_inf

        elif velocity_mode == 'file':
            velocity_file = velocity_kwargs.get('velocity_file')
            if velocity_file is None:
                raise ValueError("velocity_file must be provided for mode='file'")
            import pandas as pd
            df = pd.read_excel(velocity_file)
            time_data = df.iloc[:, 0].values
            velocity_data = df.iloc[:, 1].values
            interpolator = interp1d(time_data, velocity_data, kind='linear', fill_value='extrapolate')
            return lambda t: float(interpolator(t))

        elif velocity_mode == 'function':
            velocity_function = velocity_kwargs.get('velocity_function')
            if velocity_function is None:
                raise ValueError("velocity_function must be provided for mode='function'")
            return velocity_function

        else:
            raise ValueError(f"Invalid velocity_mode: {velocity_mode}")

    def _setup_flow_angle_profile(self, flow_angle_mode, **flow_angle_kwargs):
        if flow_angle_mode == 'constant':
            flow_angle_metocean = flow_angle_kwargs.get('flow_angle_metocean', self.flow_angle_metocean)
            flow_angle_internal = 270.0 - flow_angle_metocean
            return lambda t: flow_angle_internal

        elif flow_angle_mode == 'file':
            flow_angle_file = flow_angle_kwargs.get('flow_angle_file')
            if flow_angle_file is None:
                raise ValueError("flow_angle_file must be provided for mode='file'")
            import pandas as pd
            df = pd.read_excel(flow_angle_file)
            time_data = df.iloc[:, 0].values
            angle_data = df.iloc[:, 1].values
            interpolator = interp1d(time_data, angle_data, kind='linear', fill_value='extrapolate')
            return lambda t: 270.0 - float(interpolator(t))

        elif flow_angle_mode == 'function':
            flow_angle_function = flow_angle_kwargs.get('flow_angle_function')
            if flow_angle_function is None:
                raise ValueError("flow_angle_function must be provided for mode='function'")
            return lambda t: 270.0 - flow_angle_function(t)

        else:
            raise ValueError(f"Invalid flow_angle_mode: {flow_angle_mode}")

    def run(self, velocity_mode, total_time, progress=True, flow_angle_mode='constant', **velocity_kwargs):
        """
        Run simulation to completion.

        Parameters
        ----------
        velocity_mode : str
            'constant', 'file', or 'function'
        total_time : float
            Total simulation time [s]
        progress : bool
            Show tqdm progress bar
        flow_angle_mode : str
            'constant', 'file', or 'function'
        **velocity_kwargs
            For 'constant': U_inf=1.0
            For 'file': velocity_file='path.xlsx'
            For 'function': velocity_function=callable
            For flow_angle 'constant': flow_angle_metocean=270.0
            For flow_angle 'file': flow_angle_file='path.xlsx'
            For flow_angle 'function': flow_angle_function=callable

        Returns
        -------
        results : pd.DataFrame
            Simulation results
        """
        get_velocity = self._setup_velocity_profile(velocity_mode, **velocity_kwargs)
        get_flow_angle = self._setup_flow_angle_profile(flow_angle_mode, **velocity_kwargs)

        all_vortices = []
        time = 0.0
        shed_counters = [0 for _ in self.cylinders]

        U_inf_0 = get_velocity(time)
        Re_0 = U_inf_0 * self.D_ref / self.nu
        St_0 = compute_strouhal_number(Re_0)
        shed_period_0 = compute_shedding_period(U_inf_0, self.D_ref, St_0)
        next_shed_times = [shed_period_0 for _ in self.cylinders]

        enable_eddy_viscosity = (Re_0 > self.turbulence_thresholds['eddy_viscosity'])
        enable_stochastic_shedding = (Re_0 > self.turbulence_thresholds['stochastic_shedding'])
        enable_core_saturation = (Re_0 > self.turbulence_thresholds['core_saturation'])

        print("Starting vortex method simulation...")
        print(f"  Kinematic viscosity = {self.nu:.6e} m²/s")
        print(f"  Initial velocity U_inf(0) = {U_inf_0:.3f} m/s")
        print(f"  Initial Reynolds Re(0) = {Re_0:.2e}")
        print(f"  Strouhal number St(0) = {St_0:.3f}")
        print(f"  Turbulence modeling:")
        print(f"    Eddy viscosity: {'enabled' if enable_eddy_viscosity else 'disabled'}")
        print(f"    Stochastic shedding: {'enabled' if enable_stochastic_shedding else 'disabled'}")
        print(f"    Core saturation: {'enabled' if enable_core_saturation else 'disabled'}")
        print(f"  Initial shedding period = {shed_period_0:.3f} s\n")

        history_data = {
            'time': [],
            'U_inf': [],
            'Re': [],
            'St': [],
            'n_vortices': [],
            'vortex_field': []
        }

        for i in range(len(self.measurement_points)):
            history_data[f'probe_{i}_ux'] = []
            history_data[f'probe_{i}_uy'] = []
            history_data[f'probe_{i}_vmag'] = []

        next_save_time = self.save_interval if self.save_interval > 0 else np.inf
        num_steps = int(total_time / self.dt)

        for step in tqdm(range(num_steps), desc="Simulating", mininterval=0.5, unit="step", disable=not progress):

            for c_idx, cyl in enumerate(self.cylinders):
                if time >= next_shed_times[c_idx]:
                    U_inf_now = get_velocity(time)
                    Re_now = U_inf_now * self.D_ref / self.nu
                    St_now = compute_strouhal_number(Re_now)
                    Gamma_mag_now = compute_circulation_magnitude(U_inf_now, self.D_ref, St_now)

                    upper = (shed_counters[c_idx] % 2 == 0)
                    new_vortex = shed_vortex_with_turbulence(
                        cyl, time, upper, Re_now, U_inf_now, self.D_ref, Gamma_mag_now,
                        self.rotation_angle, get_flow_angle(time), self.theta_sep,
                        enable_stochastic=enable_stochastic_shedding
                    )
                    all_vortices.append(new_vortex)
                    shed_counters[c_idx] += 1

                    shed_period_now = compute_shedding_period(U_inf_now, self.D_ref, St_now)
                    next_shed_times[c_idx] += shed_period_now / 2  # Two vortices per Strouhal period

            N = len(all_vortices)
            if N > 0:
                birth_times = np.array([v['birth_t'] for v in all_vortices])
                ages = time - birth_times

                if enable_eddy_viscosity:
                    U_inf_now = get_velocity(time)
                    Re_now = U_inf_now * self.D_ref / self.nu

                    for i in range(N):
                        x_vortex = all_vortices[i]['x']
                        y_vortex = all_vortices[i]['y']
                        nu_eff = compute_eddy_viscosity_field(x_vortex, y_vortex, self.D_ref, self.nu, U_inf_now, Re_now)
                        age = ages[i]
                        sigma_0_vortex = all_vortices[i]['sigma_0']
                        sigma_new = np.sqrt(sigma_0_vortex**2 + 4 * nu_eff * age)

                        if enable_core_saturation:
                            sigma_max = self.sigma_max_factor * self.D_ref
                            sigma_new = min(sigma_new, sigma_max)

                        all_vortices[i]['sigma'] = sigma_new
                else:
                    for i in range(N):
                        sigma_0_vortex = all_vortices[i]['sigma_0']
                        age = ages[i]
                        sigma_new = np.sqrt(sigma_0_vortex**2 + 4 * self.nu * age)
                        all_vortices[i]['sigma'] = sigma_new

                for i in range(N):
                    t_age = time - all_vortices[i]['birth_t']
                    growth_factor = circulation_growth_sigmoid(t_age, all_vortices[i]['T_form'])
                    all_vortices[i]['gamma'] = all_vortices[i]['gamma_target'] * growth_factor

                sources_x = np.array([v['x'] for v in all_vortices])
                sources_y = np.array([v['y'] for v in all_vortices])
                gammas = np.array([v['gamma'] for v in all_vortices])
                sigmas = np.array([v['sigma'] for v in all_vortices])

                U_inf_t = get_velocity(time)
                U1, V1 = compute_velocities_at_points(sources_x, sources_y, sources_x, sources_y, gammas, sigmas, U_inf_t, get_flow_angle(time))
                k1x = self.dt * U1
                k1y = self.dt * V1

                U_inf_mid = get_velocity(time + 0.5*self.dt)
                mid_x = sources_x + 0.5 * k1x
                mid_y = sources_y + 0.5 * k1y
                U2, V2 = compute_velocities_at_points(mid_x, mid_y, sources_x, sources_y, gammas, sigmas, U_inf_mid, get_flow_angle(time + 0.5*self.dt))
                k2x = self.dt * U2
                k2y = self.dt * V2

                mid_x = sources_x + 0.5 * k2x
                mid_y = sources_y + 0.5 * k2y
                U3, V3 = compute_velocities_at_points(mid_x, mid_y, sources_x, sources_y, gammas, sigmas, U_inf_mid, get_flow_angle(time + 0.5*self.dt))
                k3x = self.dt * U3
                k3y = self.dt * V3

                U_inf_end = get_velocity(time + self.dt)
                end_x = sources_x + k3x
                end_y = sources_y + k3y
                U4, V4 = compute_velocities_at_points(end_x, end_y, sources_x, sources_y, gammas, sigmas, U_inf_end, get_flow_angle(time + self.dt))
                k4x = self.dt * U4
                k4y = self.dt * V4

                new_x = sources_x + (k1x + 2*k2x + 2*k3x + k4x) / 6
                new_y = sources_y + (k1y + 2*k2y + 2*k3y + k4y) / 6
                for i in range(N):
                    all_vortices[i]['x'] = new_x[i]
                    all_vortices[i]['y'] = new_y[i]

                all_vortices = [v for v in all_vortices if v['x'] <= self.x_removal]

            U_inf_now = get_velocity(time)
            Re_now = U_inf_now * self.D_ref / self.nu
            St_now = compute_strouhal_number(Re_now)

            history_data['time'].append(time)
            history_data['U_inf'].append(U_inf_now)
            history_data['Re'].append(Re_now)
            history_data['St'].append(St_now)
            history_data['n_vortices'].append(len(all_vortices))

            if len(all_vortices) > 0:
                sources_x = np.array([v['x'] for v in all_vortices])
                sources_y = np.array([v['y'] for v in all_vortices])
                gammas = np.array([v['gamma'] for v in all_vortices])
                sigmas = np.array([v['sigma'] for v in all_vortices])

                if len(self.measurement_points) > 0:
                    probe_x = np.array([p[0] for p in self.measurement_points])
                    probe_y = np.array([p[1] for p in self.measurement_points])
                    ux, uy = compute_velocities_at_points(
                        probe_x, probe_y, sources_x, sources_y, gammas, sigmas, U_inf_now, get_flow_angle(time)
                    )
                    for i in range(len(self.measurement_points)):
                        history_data[f'probe_{i}_ux'].append(ux[i])
                        history_data[f'probe_{i}_uy'].append(uy[i])
                        history_data[f'probe_{i}_vmag'].append(np.sqrt(ux[i]**2 + uy[i]**2))
                else:
                    for i in range(len(self.measurement_points)):
                        history_data[f'probe_{i}_ux'].append(U_inf_now * np.cos(np.radians(get_flow_angle(time))))
                        history_data[f'probe_{i}_uy'].append(U_inf_now * np.sin(np.radians(get_flow_angle(time))))
                        history_data[f'probe_{i}_vmag'].append(U_inf_now)
            else:
                for i in range(len(self.measurement_points)):
                    history_data[f'probe_{i}_ux'].append(U_inf_now * np.cos(np.radians(get_flow_angle(time))))
                    history_data[f'probe_{i}_uy'].append(U_inf_now * np.sin(np.radians(get_flow_angle(time))))
                    history_data[f'probe_{i}_vmag'].append(U_inf_now)

            if time >= next_save_time or (self.save_interval == 0 and step == num_steps - 1):
                vortex_snapshot = [
                    {
                        'x': v['x'],
                        'y': v['y'],
                        'gamma': v['gamma'],
                        'sigma': v['sigma'],
                        'birth_t': v['birth_t'],
                        'age': time - v['birth_t']
                    }
                    for v in all_vortices
                ]
                history_data['vortex_field'].append(vortex_snapshot)
                next_save_time += self.save_interval
            else:
                history_data['vortex_field'].append(None)

            time += self.dt

        print(f"Simulation complete. Total vortices in domain: {len(all_vortices)}")

        self.results_df = pd.DataFrame(history_data)
        return self.results_df

    def save_results(self, filename):
        """Save results DataFrame to pickle file"""
        if self.results_df is None:
            raise ValueError("No results to save. Run simulation first.")
        self.results_df.to_pickle(filename)
        print(f"Results saved to {filename}")

    @staticmethod
    def load_results(filename):
        """Load results DataFrame from pickle file"""
        return pd.read_pickle(filename)
