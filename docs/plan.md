# Dynamic Vortex Tracking Implementation Plan

## Core Changes Required

Your script currently uses **prescribed spatial vortex positions** via formulas. You need to replace this with **time-integrated vortex trajectories** where vortices are shed periodically and advected by the velocity field.

---

## 1. Replace Static Vortex Data Structure

**Current:** `get_vortex_positions()` returns positions from formulas

**New:** Maintain a dynamic list of `Vortex` objects

```python
@dataclass
class Vortex:
    x: float           # Current x position [D]
    y: float           # Current y position [D]  
    gamma: float       # Circulation (+ or -)
    t_birth: float     # Time when shed [T]
    
    def sigma(self, t, params):
        """Core radius via Lamb-Oseen: σ²(t) = σ₀² + 4ν(t - t_birth)"""
        age = t - self.t_birth
        nu = params['U_inf'] * params['D'] / params['Re']
        return np.sqrt(params['sigma_0']**2 + 4 * nu * age)
```

---

## 2. Vortex Shedding Parameters

Add to `params` dictionary:

```python
'St': 0.2,              # Strouhal number (shedding frequency)
'x_shed': 1.0,          # Shedding x-location [D] (behind cylinder)
'y_shed': 0.5,          # Shedding y-offset [D] (half-width)
't_max': 100.0,         # Simulation time [T = D/U_inf]
'dt': 0.02,             # Integration timestep [T]
```

**Circulation per vortex (empirical):**
$$\Gamma_{\text{shed}} = 2\pi \cdot \text{St} \cdot U_\infty \cdot D$$

This follows from the definition $\text{St} = f D / U_\infty$ and the relation $\Gamma \sim U_\infty D$ per shedding cycle.

---

## 3. Vortex Shedding Logic

**Replace `get_vortex_positions()` calls with:**

```python
def shed_vortex(vortices, t, params):
    """Shed vortex if period elapsed; alternate sign"""
    T_shed = 1.0 / (params['St'] * params['U_inf'] / params['D'])
    shed_count = int(t / T_shed)
    
    # Check if new vortex due
    if t >= shed_count * T_shed and not hasattr(shed_vortex, f'shed_{shed_count}'):
        setattr(shed_vortex, f'shed_{shed_count}', True)  # Mark as shed
        
        # Alternate upper/lower
        y = params['y_shed'] if shed_count % 2 == 0 else -params['y_shed']
        gamma_sign = +1 if shed_count % 2 == 0 else -1
        
        gamma_shed = 2 * np.pi * params['St'] * params['U_inf'] * params['D']
        
        vortex = Vortex(
            x=params['x_shed'],
            y=y,
            gamma=gamma_sign * gamma_shed,
            t_birth=t
        )
        vortices.append(vortex)
```

---

## 4. Velocity Field from Dynamic Vortices

**Modify `gaussian_velocity_contribution()` to accept time-dependent σ:**

```python
def gaussian_velocity_contribution_dynamic(X, Y, vortex, t, params):
    """Compute velocity from single vortex with time-dependent core"""
    sigma = vortex.sigma(t, params)  # Core grows with age
    dx = X - vortex.x
    dy = Y - vortex.y
    r_squared = dx**2 + dy**2
    r = np.sqrt(r_squared + 1e-10)
    
    exp_factor = np.exp(-r_squared / sigma**2)
    velocity_magnitude = (vortex.gamma / (2 * np.pi * r)) * (1 - exp_factor)
    
    U = velocity_magnitude * (-dy / r)
    V = velocity_magnitude * (dx / r)
    
    return U, V

def compute_velocity_field_dynamic(X, Y, vortices, t, params):
    """Total velocity from all vortices at time t"""
    U = np.full_like(X, params['U_inf'])
    V = np.zeros_like(X)
    
    for vortex in vortices:
        U_v, V_v = gaussian_velocity_contribution_dynamic(X, Y, vortex, t, params)
        U += U_v
        V += V_v
    
    return U, V
```

---

## 5. Vortex Advection (The Integration)

**Euler method (first-order, simple):**

$$\mathbf{x}_j^{n+1} = \mathbf{x}_j^n + \mathbf{u}(\mathbf{x}_j^n, t^n) \cdot \Delta t$$

```python
def advect_vortices_euler(vortices, t, dt, params, epsilon=0.01):
    """
    Integrate vortex motion: dx/dt = u(x,t)
    
    Args:
        epsilon: Self-induction cutoff radius [D] to avoid singularities
    """
    # Store current positions (to avoid using updated positions in same step)
    positions = [(v.x, v.y) for v in vortices]
    
    for i, vortex in enumerate(vortices):
        # Compute velocity AT this vortex's position
        u = params['U_inf']
        v = 0.0
        
        for j, other in enumerate(vortices):
            if i == j:
                continue  # Skip self-induction
            
            dx = positions[i][0] - other.x
            dy = positions[i][1] - other.y
            r_sq = dx**2 + dy**2
            
            # Cutoff to prevent singularity
            r_cutoff = epsilon * params['D']
            if r_sq < r_cutoff**2:
                continue
            
            # Lamb-Oseen kernel at vortex position
            r = np.sqrt(r_sq)
            sigma_j = other.sigma(t, params)
            exp_factor = np.exp(-r_sq / sigma_j**2)
            vel_mag = (other.gamma / (2 * np.pi * r)) * (1 - exp_factor)
            
            u += vel_mag * (-dy / r)
            v += vel_mag * (dx / r)
        
        # Euler step
        vortex.x += u * dt
        vortex.y += v * dt
```

**Runge-Kutta 4th order (more accurate, recommended):**

$$\mathbf{x}^{n+1} = \mathbf{x}^n + \frac{\Delta t}{6}(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4)$$

where:
$$\mathbf{k}_1 = \mathbf{u}(t^n, \mathbf{x}^n)$$
$$\mathbf{k}_2 = \mathbf{u}(t^n + \Delta t/2, \mathbf{x}^n + \mathbf{k}_1 \Delta t/2)$$
$$\mathbf{k}_3 = \mathbf{u}(t^n + \Delta t/2, \mathbf{x}^n + \mathbf{k}_2 \Delta t/2)$$
$$\mathbf{k}_4 = \mathbf{u}(t^n + \Delta t, \mathbf{x}^n + \mathbf{k}_3 \Delta t)$$

```python
def compute_vortex_velocity(vortex, vortices, t, params, epsilon=0.01):
    """Helper: velocity at single vortex position"""
    u = params['U_inf']
    v = 0.0
    
    for other in vortices:
        if vortex is other:
            continue
        
        dx = vortex.x - other.x
        dy = vortex.y - other.y
        r_sq = dx**2 + dy**2
        
        r_cutoff = epsilon * params['D']
        if r_sq < r_cutoff**2:
            continue
        
        r = np.sqrt(r_sq)
        sigma = other.sigma(t, params)
        exp_factor = np.exp(-r_sq / sigma**2)
        vel_mag = (other.gamma / (2 * np.pi * r)) * (1 - exp_factor)
        
        u += vel_mag * (-dy / r)
        v += vel_mag * (dx / r)
    
    return np.array([u, v])

def advect_vortices_rk4(vortices, t, dt, params):
    """4th-order Runge-Kutta integration"""
    x0 = np.array([[v.x, v.y] for v in vortices])
    
    # Stage 1
    k1 = np.array([compute_vortex_velocity(v, vortices, t, params) for v in vortices])
    
    # Stage 2 (temporary positions)
    for i, v in enumerate(vortices):
        v.x, v.y = x0[i] + 0.5 * dt * k1[i]
    k2 = np.array([compute_vortex_velocity(v, vortices, t + 0.5*dt, params) for v in vortices])
    
    # Stage 3
    for i, v in enumerate(vortices):
        v.x, v.y = x0[i] + 0.5 * dt * k2[i]
    k3 = np.array([compute_vortex_velocity(v, vortices, t + 0.5*dt, params) for v in vortices])
    
    # Stage 4
    for i, v in enumerate(vortices):
        v.x, v.y = x0[i] + dt * k3[i]
    k4 = np.array([compute_vortex_velocity(v, vortices, t + dt, params) for v in vortices])
    
    # Final update
    for i, v in enumerate(vortices):
        pos_new = x0[i] + (dt / 6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
        v.x, v.y = pos_new
```

---

## 6. Main Time-Integration Loop

**Replace `compute_time_varying_velocity()` with:**

```python
def simulate_dynamic_vortex_street(params):
    """
    Main simulation loop: shed vortices, advect them, sample velocity
    
    Returns:
        times: Time array [T]
        vortices_history: List of vortex lists at each timestep
        u_probe, v_probe: Velocity timeseries at probe location
    """
    dt = params['dt']
    t_max = params['t_max']
    n_steps = int(t_max / dt)
    
    # Probe location
    x_probe = params['probe_x']
    y_probe = params['probe_y']
    
    # Storage
    times = []
    vortices_history = []
    u_probe_history = []
    v_probe_history = []
    
    # Initialize
    vortices = []
    t = 0.0
    
    # Time loop
    for step in range(n_steps):
        # Shed vortex if needed
        shed_vortex(vortices, t, params)
        
        # Advect vortices (choose method)
        # advect_vortices_euler(vortices, t, dt, params)
        advect_vortices_rk4(vortices, t, dt, params)
        
        # Remove vortices outside domain
        vortices = [v for v in vortices 
                   if v.x < params['x_max'] and abs(v.y) < params['y_max']]
        
        # Sample velocity at probe
        u_p = params['U_inf']
        v_p = 0.0
        for vortex in vortices:
            dx = x_probe - vortex.x
            dy = y_probe - vortex.y
            r_sq = dx**2 + dy**2 + 1e-10
            r = np.sqrt(r_sq)
            sigma = vortex.sigma(t, params)
            exp_factor = np.exp(-r_sq / sigma**2)
            vel_mag = (vortex.gamma / (2 * np.pi * r)) * (1 - exp_factor)
            u_p += vel_mag * (-dy / r)
            v_p += vel_mag * (dx / r)
        
        # Store
        times.append(t)
        vortices_history.append([Vortex(v.x, v.y, v.gamma, v.t_birth) for v in vortices])
        u_probe_history.append(u_p)
        v_probe_history.append(v_p)
        
        # Advance time
        t += dt
        
        if step % 100 == 0:
            print(f"Step {step}/{n_steps}, t={t:.2f}, N_vortices={len(vortices)}")
    
    return (np.array(times), vortices_history, 
            np.array(u_probe_history), np.array(v_probe_history))
```

---

## 7. Update Visualization Functions

**Modify Plot 4 (time-varying velocity):**

Replace the analytical translation model with real simulation results:

```python
# In create_visualization(), replace Plot 4 section with:

times, vortices_history, u_probe, v_probe = simulate_dynamic_vortex_street(params)
vel_mag_probe = np.sqrt(u_probe**2 + v_probe**2)

ax4.plot(times, u_probe, 'b-', label='u (streamwise)', linewidth=2)
ax4.plot(times, v_probe, 'r-', label='v (transverse)', linewidth=2)
ax4.plot(times, vel_mag_probe, 'g--', label='|v|', linewidth=1.5, alpha=0.7)

# Compute actual shedding period from data (FFT or peak detection)
from scipy import signal
peaks, _ = signal.find_peaks(u_probe, height=0)
if len(peaks) > 1:
    T_measured = np.mean(np.diff(times[peaks]))
    ax4.axvline(T_measured, color='black', linestyle=':', 
                label=f'Measured period T={T_measured:.2f}')
```

**Add animation of vortex convection:**

```python
from matplotlib.animation import FuncAnimation

def animate_vortex_street(vortices_history, times, params):
    """Animate vortices moving downstream"""
    fig, ax = plt.subplots(figsize=(14, 6))
    
    def update(frame):
        ax.clear()
        vortices = vortices_history[frame]
        t = times[frame]
        
        # Separate by sign
        x_pos = [v.x for v in vortices if v.gamma > 0]
        y_pos = [v.y for v in vortices if v.gamma > 0]
        x_neg = [v.x for v in vortices if v.gamma < 0]
        y_neg = [v.y for v in vortices if v.gamma < 0]
        
        ax.scatter(x_pos, y_pos, c='red', marker='o', s=50, label='+Γ', zorder=5)
        ax.scatter(x_neg, y_neg, c='blue', marker='o', s=50, label='-Γ', zorder=5)
        
        # Draw cores
        for v in vortices[::3]:  # Every 3rd to reduce clutter
            sigma = v.sigma(t, params)
            circle = plt.Circle((v.x, v.y), sigma, fill=False, 
                              edgecolor='gray', linestyle=':', alpha=0.3)
            ax.add_patch(circle)
        
        # Cylinder
        cyl = plt.Circle((0, 0), 0.5*params['D'], color='gray', zorder=10)
        ax.add_patch(cyl)
        
        ax.set_xlim([-1, 15])
        ax.set_ylim([-3, 3])
        ax.set_xlabel('x/D')
        ax.set_ylabel('y/D')
        ax.set_title(f'Vortex Street Evolution - t = {t:.2f}')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    anim = FuncAnimation(fig, update, frames=len(times), interval=50)
    return anim
```

---

## 8. Critical Numerical Parameters

| Parameter | Typical Value | Constraint | Reason |
|-----------|---------------|------------|--------|
| **dt** | 0.01 - 0.05 | dt < σ₀/U_inf | CFL stability |
| **epsilon** | 0.01D | δ << σ₀ | Self-induction cutoff |
| **x_shed** | 0.5 - 1.5D | Wake formation length | Empirical from CFD |
| **y_shed** | 0.4 - 0.6D | ~ D/2 | Initial separation |
| **Γ_shed** | 2πSt·U·D | St ≈ 0.2 | Kármán relation |

**Timestep criterion (CFL-like):**
$$\Delta t < \frac{\sigma_0}{|\mathbf{u}|_{\max}} \approx \frac{0.1D}{U_\infty}$$

For D = 1, U = 1: dt < 0.1, use dt ≈ 0.02 for safety.

---

## 9. Physics Validation Checks

After running simulation, verify:

1. **Strouhal number from data:**
$$\text{St}_{\text{measured}} = f_{\text{peak}} \cdot D / U_\infty \approx 0.2$$
(Use FFT of u_probe to find peak frequency)

2. **Wake spreading matches model:**
Extract y-positions of vortices vs x, fit to:
$$h(x) = h_0(1 + \alpha x/D)$$
Compare α with `spreading_coef` parameter.

3. **Core growth:**
Measure σ from vorticity peaks, verify:
$$\sigma^2(t) = \sigma_0^2 + 4\nu t$$

---

## 10. Implementation Sequence

1. **Add Vortex dataclass** to top of script
2. **Add shedding parameters** to `params` dictionary
3. **Create `shed_vortex()`** function
4. **Modify velocity functions** to use vortex list + time
5. **Implement `advect_vortices_rk4()`** (start with Euler if debugging)
6. **Create `simulate_dynamic_vortex_street()`** main loop
7. **Update Plot 4** to use real simulation data
8. **Add animation** (optional but impressive)
9. **Run validation checks** on output

**Expected outcome:** Oscillating velocity signal at probe, vortices convecting downstream, growing cores, diverging wake - all from **first-principles integration**, not prescribed geometry.