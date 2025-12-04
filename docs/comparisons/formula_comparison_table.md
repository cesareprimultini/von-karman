# FORMULA-BY-FORMULA COMPARISON TABLE
## Von Kármán Vortex Street Simulations

---

## 1. VORTEX CORE RADIUS σ(t) - Viscous Diffusion

### Mathematical Formula:
```
σ²(t) = σ₀² + 4ν·t_age

where:
  σ₀ = initial core radius
  ν = kinematic viscosity = U∞·D/Re
  t_age = age of vortex (time since birth or x/U∞)
```

### Implementation Comparison:

| Script | Age Definition | Implementation | Comments |
|--------|----------------|----------------|----------|
| **von_karman** | `age = t - t_birth` | `sigma_sq = sigma_0**2 + 4*nu*age` | TRUE temporal age |
| **simple_vortex** | `age = t - t_birth` | `sigma = sqrt(sigma_0**2 + 4*nu*age)` | TRUE temporal age |
| **test.py** | `age = (x - x_cyl)/U_inf` | `age_up = (x_up - x_cyl)/U_inf` | SPATIAL age (kinematic) |
| **test01.py** | `age = t - birth_t` | `sigma = sqrt(sigma_0**2 + 4*nu*age)` | TRUE temporal age |

### Key Difference:
- **test.py uses SPATIAL age**: Assumes vortices convect at U∞, so age = distance/velocity
- **Others use TEMPORAL age**: Track actual time since vortex was shed
- **Physics**: Temporal age is correct; spatial is approximation valid only if U_vortex ≈ U∞

---

## 2. LAMB-OSEEN VELOCITY FIELD

### Mathematical Formula:
```
Induced velocity at point (x,y) from vortex at (x_v, y_v):

r = √[(x - x_v)² + (y - y_v)²]
V_θ(r) = (Γ/2πr) × [1 - exp(-r²/σ²)]

Cartesian components:
u = -V_θ × (y - y_v)/r  [x-component]
v = +V_θ × (x - x_v)/r  [y-component]
```

### Implementation Comparison:

#### **von_karman_simulator.py:**
```python
dx = X - x_v
dy = Y - y_v
r_sq = dx**2 + dy**2
r = np.sqrt(r_sq + 1e-10)

V_theta = (gamma / (2*np.pi*r)) * (1 - np.exp(-r_sq/(sigma**2 + 1e-10)))

U = -V_theta * dy / r
V = V_theta * dx / r
```

#### **simple_vortex_street.py:**
```python
dx = x - x_vortex
dy = y - y_vortex
r_sq = dx**2 + dy**2
r = np.sqrt(r_sq + 1e-10)

exp_factor = np.exp(-r_sq / (sigma**2 + 1e-10))
vel_mag = (gamma / (2*np.pi*r)) * (1 - exp_factor)

u = -vel_mag * dy / r
v = vel_mag * dx / r
```

#### **test.py:**
```python
dx = x - vortex['x']
dy = y - vortex['y']
r = np.sqrt(dx**2 + dy**2)

if r < 1e-3:
    continue

sigma = vortex['sigma']
exp_factor = np.exp(-r**2 / sigma**2)
V_theta = (vortex['gamma'] / (2*np.pi*r)) * (1 - exp_factor)

u += -V_theta * dy / r
v += +V_theta * dx / r
```

#### **test01.py:**
```python
def induced_velocity(dx, dy, r, gamma, sigma):
    if r < 1e-6:
        return 0.0, 0.0
    
    exp_factor = np.exp(-r**2 / sigma**2)
    V_theta = (gamma / (2*np.pi*r)) * (1 - exp_factor)
    
    u_ind = -V_theta * dy / r
    v_ind = V_theta * dx / r
    return u_ind, v_ind
```

### Differences:
- **Numerical stability**: All use small epsilon (1e-10, 1e-6, 1e-3)
- **Cutoff methods**: Distance check vs epsilon addition
- **Result**: MATHEMATICALLY IDENTICAL

---

## 3. CIRCULATION MAGNITUDE Γ

### Mathematical Formula:
```
Strouhal-based:
  Γ = C_γ × St × U∞ × D

where:
  St = Strouhal number ≈ 0.2 for Re ~ 100-200
  C_γ = circulation coefficient (typically 2π for dimensional consistency)
  U∞ = freestream velocity
  D = cylinder diameter
```

### Implementation Comparison:

| Script | Formula | C_γ | St | Value |
|--------|---------|-----|----|----|
| **von_karman** | `C_gamma * St * U_inf * D` | 2π | 0.2 | ~1.26·U∞·D |
| **simple_vortex** | `2.5 * U_inf * D` | 2.5 | N/A | 2.5·U∞·D |
| **test.py** | `2*pi * St * U_inf * D` | 2π | 0.2 | ~1.26·U∞·D |
| **test01.py** | `2*pi * St * U_inf * D` | 2π | 0.2 | ~1.26·U∞·D |

### Code Snippets:

**von_karman:**
```python
C_gamma = 2 * np.pi  # In params
Gamma = params['C_gamma'] * params['St'] * params['U_inf'] * params['D']
```

**simple_vortex:**
```python
Gamma_shed = 2.5 * U_inf * D  # Empirical constant
```

**test.py:**
```python
St = 0.2
Gamma = 2.0 * np.pi * St * U_inf * D
```

**test01.py:**
```python
St = 0.2
Gamma_mag = 2.0 * np.pi * St * U_inf * D
```

### Analysis:
- **simple_vortex**: Uses empirical 2.5 (not physics-based, about 2× too high)
- **Others**: Use standard 2π·St ≈ 1.26 (correct dimensionally)
- **Impact**: Higher Γ → stronger wake, faster vortex rotation

---

## 4. TIME INTEGRATION / VORTEX ADVECTION

### A. EULER METHOD (First-Order)

#### Mathematical Formula:
```
x_{n+1} = x_n + v_x(x_n, y_n)·dt
y_{n+1} = y_n + v_y(x_n, y_n)·dt

Error per step: O(dt²)
Total error: O(dt)
```

#### **von_karman_simulator.py:**
```python
def advect_vortices_euler(vortices, params, t, dt):
    new_vortices = []
    for vortex in vortices:
        # Compute velocity at current position
        u, v = compute_velocity_at_point(vortex.x, vortex.y, vortices, t, params)
        
        # Update position
        x_new = vortex.x + u * dt
        y_new = vortex.y + v * dt
        
        new_vortices.append(Vortex(x_new, y_new, vortex.gamma, vortex.t_birth))
    return new_vortices
```

#### **simple_vortex_street.py:**
```python
new_vortices = []
for xv, yv, gamma, t_birth in vortices:
    u, v = get_velocity_at_point(xv, yv, vortices, t)
    
    xv_new = xv + u * dt
    yv_new = yv + v * dt
    
    # Domain check
    if xv_new < 15*D and abs(yv_new) < 4*D:
        new_vortices.append([xv_new, yv_new, gamma, t_birth])

vortices = new_vortices
```

### B. RK4 METHOD (Fourth-Order)

#### Mathematical Formula:
```
k₁ = v(x_n) · dt
k₂ = v(x_n + k₁/2) · dt
k₃ = v(x_n + k₂/2) · dt
k₄ = v(x_n + k₃) · dt

x_{n+1} = x_n + (k₁ + 2k₂ + 2k₃ + k₄)/6

Error per step: O(dt⁵)
Total error: O(dt⁴)
```

#### **test01.py:**
```python
def rk4_advect(vortex_idx, vortices, dt):
    def vel_func(pos):
        return velocity_at_point(pos[0], pos[1], vortices, 
                                include_images=True, 
                                exclude_self_idx=vortex_idx)
    
    x0, y0 = vortices[vortex_idx]['x'], vortices[vortex_idx]['y']
    pos0 = np.array([x0, y0])
    
    # Four stages of RK4
    k1 = np.array(vel_func(pos0)) * dt
    k2 = np.array(vel_func(pos0 + 0.5*k1)) * dt
    k3 = np.array(vel_func(pos0 + 0.5*k2)) * dt
    k4 = np.array(vel_func(pos0 + k3)) * dt
    
    # Weighted average
    pos_new = pos0 + (k1 + 2*k2 + 2*k3 + k4) / 6
    
    return pos_new[0], pos_new[1]

# In main loop:
new_positions = []
for idx in range(len(all_vortices)):
    new_x, new_y = rk4_advect(idx, all_vortices, dt)
    new_positions.append((new_x, new_y))

for idx, (new_x, new_y) in enumerate(new_positions):
    all_vortices[idx]['x'] = new_x
    all_vortices[idx]['y'] = new_y
```

### C. NO ADVECTION (Algebraic)

#### **test.py:**
```python
# Vortices at fixed positions
x_up = x_start + k * l
y_up = y_cyl + h/2

# Assumes convection at U∞ with no dynamics
```

### Accuracy Comparison:

| Method | Error/Step | Total Error | Computational Cost | Scripts |
|--------|-----------|-------------|-------------------|---------|
| **None** | N/A | N/A (static) | O(1) | test.py |
| **Euler** | O(dt²) | O(dt) | O(N²) per step | von_karman, simple_vortex |
| **RK4** | O(dt⁵) | O(dt⁴) | O(4N²) per step | test01.py |

### Stability Requirements:

**Euler Method:**
```
dt < σ₀/U_max ≈ 0.1·D/U∞
```
For D=1, U∞=1, σ₀=0.1 → **dt < 0.1**

**RK4 Method:**
```
dt < 2·σ₀/U_max ≈ 0.2·D/U∞
```
Can use **2× larger timestep** for same stability

---

## 5. IMAGE VORTICES (Boundary Conditions)

### Mathematical Theory:

For a cylinder of radius `a` centered at (x_c, y_c):
```
Complex notation:
  z_img = z_c + a²/(z - z_c)*

Cartesian coordinates:
  r_rel² = (x - x_c)² + (y - y_c)²
  
  x_img = x_c + a²·(x - x_c)/r_rel²
  y_img = y_c + a²·(y - y_c)/r_rel²
  
  Γ_img = -Γ  (opposite sign!)
  σ_img = σ   (same core radius)
```

### Physics:
- Image vortex enforces **no-penetration** boundary condition
- Ensures **v·n = 0** on cylinder surface
- Required for **physical accuracy** near boundaries

### Implementation Comparison:

#### **von_karman_simulator.py:**
```python
def get_image_vortex_position(vortex, params):
    """
    Single cylinder at origin.
    """
    x, y = vortex.x, vortex.y
    a = params['D'] / 2  # Cylinder radius
    
    r_sq = x**2 + y**2
    
    # Validity check
    if r_sq <= a**2:
        return None  # Vortex inside cylinder (invalid)
    
    x_img = (a**2 * x) / r_sq
    y_img = (a**2 * y) / r_sq
    
    return (x_img, y_img)

# Usage in velocity calculation:
if params.get('use_image_vortices', False):
    image_pos = get_image_vortex_position(vortex, params)
    if image_pos is not None:
        x_img, y_img = image_pos
        sigma_img = vortex.sigma(t, params)
        gamma_img = -vortex.gamma  # OPPOSITE SIGN
        
        U_img, V_img = gaussian_velocity_contribution(
            X, Y, x_img, y_img, sigma_img, gamma_img
        )
        U += U_img
        V += V_img
```

#### **test01.py:**
```python
def get_image_vortex(vortex, cyl):
    """
    Multiple cylinders, relative coordinates.
    """
    # Relative position
    x_rel = vortex['x'] - cyl['x']
    y_rel = vortex['y'] - cyl['y']
    r_sq = x_rel**2 + y_rel**2
    
    # Validity check
    if r_sq <= a**2 + 1e-6:
        return None
    
    # Image in relative coordinates
    x_img_rel = (a**2 * x_rel) / r_sq
    y_img_rel = (a**2 * y_rel) / r_sq
    
    # Convert to absolute coordinates
    x_img = cyl['x'] + x_img_rel
    y_img = cyl['y'] + y_img_rel
    
    gamma_img = -vortex['gamma']  # OPPOSITE SIGN
    sigma_img = vortex['sigma']
    
    return {'x': x_img, 'y': y_img, 'gamma': gamma_img, 'sigma': sigma_img}

# Usage in velocity calculation:
for cyl in cylinders:
    img = get_image_vortex(vortex, cyl)
    if img:
        dx_img = x - img['x']
        dy_img = y - img['y']
        r_img = np.sqrt(dx_img**2 + dy_img**2)
        u_img, v_img = induced_velocity(dx_img, dy_img, r_img, 
                                       img['gamma'], img['sigma'])
        u += u_img
        v += v_img
```

#### **simple_vortex_street.py & test.py:**
```python
# NO IMAGE VORTICES IMPLEMENTED
# Boundary condition violated!
```

### Key Differences:

| Feature | von_karman | test01.py | simple_vortex | test.py |
|---------|-----------|-----------|---------------|---------|
| **Implemented?** | YES (optional) | YES (always) | NO | NO |
| **Cylinders** | 1 (at origin) | 4 (grid) | N/A | N/A |
| **Coordinates** | Absolute | Relative | N/A | N/A |
| **Images per vortex** | 1 | 4 | 0 | 0 |
| **BC accuracy** | HIGH | HIGH | NONE | NONE |

### Computational Cost:

**Without images:**
```
N vortices → N² velocity evaluations
```

**With images (single cylinder):**
```
N vortices → N² + N = N(N+1) evaluations
```

**With images (4 cylinders):**
```
N vortices → N² + 4N = N(N+4) evaluations
```

For N=100: 10,000 → 10,400 (4% increase)

---

## 6. VORTEX SHEDDING

### Mathematical Formula:

**Shedding Period:**
```
T_shed = D / (St × U∞)

where St ≈ 0.2 for Re ~ 100-200
```

**Shedding Condition:**
```
t ≥ n·T_shed → shed vortex n
```

**Separation Angle:**
```
θ_sep ≈ 70-80° from rear stagnation point
θ_upper = π - θ_sep
θ_lower = π + θ_sep
```

**Shed Position:**
```
x_shed = (D/2 + σ₀)·cos(θ)
y_shed = (D/2 + σ₀)·sin(θ)
```

**Circulation Sign:**
```
Upper row: Γ_upper = -|Γ|  (clockwise, negative)
Lower row: Γ_lower = +|Γ|  (counter-clockwise, positive)
```

### Implementation Comparison:

#### **Dynamic Shedding (von_karman, simple_vortex, test01):**

**von_karman & test01:**
```python
T_shed = D / (St * U_inf)
shed_period = D / (St * U_inf)

# In time loop:
for c_idx, cyl in enumerate(cylinders):
    if time >= next_shed_times[c_idx]:
        # Alternate upper/lower
        upper = (shed_counters[c_idx] % 2 == 0)
        
        # Separation angle
        sign = 1 if upper else -1
        theta = np.pi - sign * theta_sep
        
        # Position just outside cylinder
        x_shed = cyl['x'] + a * np.cos(theta)
        y_shed = cyl['y'] + a * np.sin(theta)
        
        # Circulation (upper is clockwise = negative)
        gamma = -Gamma_mag if upper else Gamma_mag
        
        new_vortex = {'x': x_shed, 'y': y_shed, 
                     'gamma': gamma, 'sigma': sigma_0, 
                     'birth_t': time}
        all_vortices.append(new_vortex)
        
        shed_counters[c_idx] += 1
        next_shed_times[c_idx] += shed_period
```

**simple_vortex_street.py:**
```python
T_shed = D / (St * U_inf)

# In time loop:
expected_count = int(t / T_shed)

while shed_count < expected_count:
    is_upper = (shed_count % 2 == 0)
    
    # Separation angle
    theta = np.pi - theta_sep if is_upper else np.pi + theta_sep
    
    # Position
    x_shed = (D/2 + sigma_0) * np.cos(theta)
    y_shed = (D/2 + sigma_0) * np.sin(theta)
    
    # Circulation
    gamma = -Gamma_shed if is_upper else Gamma_shed
    
    vortices.append([x_shed, y_shed, gamma, t])
    shed_count += 1
```

#### **Algebraic Positioning (test.py):**

```python
# Pre-position ALL vortices at once
h = 0.281 * D  # von Kármán stability
l = 1.0 * D    # Wavelength
a = 0.5 * l    # Stagger

for k in range(n_pairs):
    # Upper row: positive circulation
    x_up = x_start + k * l
    y_up = y_cyl + h/2
    age_up = (x_up - x_cyl) / U_inf
    sigma_up = sqrt(sigma_0**2 + 4*nu*age_up)
    vortices.append({'x': x_up, 'y': y_up, 
                    'gamma': +Gamma, 'sigma': sigma_up})
    
    # Lower row: negative circulation
    x_low = x_start + a + k * l
    y_low = y_cyl - h/2
    age_low = (x_low - x_cyl) / U_inf
    sigma_low = sqrt(sigma_0**2 + 4*nu*age_low)
    vortices.append({'x': x_low, 'y': y_low, 
                    'gamma': -Gamma, 'sigma': sigma_low})
```

### Shedding Parameters:

| Script | Method | Separation Angle | Sign Convention | Multiple Cylinders |
|--------|--------|------------------|-----------------|-------------------|
| **von_karman** | Dynamic | θ_sep = 1.4 rad (80°) | Upper: -, Lower: + | No |
| **simple_vortex** | Dynamic | θ_sep = 70° | Upper: -, Lower: + | No |
| **test.py** | Algebraic | N/A | Upper: +, Lower: - | YES (4) |
| **test01.py** | Dynamic | θ_sep = 1.4 rad (80°) | Upper: -, Lower: + | YES (4) |

### Key Differences:

1. **Dynamic vs Static:**
   - Dynamic: Vortices shed one at a time based on Strouhal number
   - Static: All vortices pre-positioned in space

2. **Sign Convention:**
   - Dynamic scripts (von_karman, simple_vortex, test01): Upper=-Γ, Lower=+Γ
   - Static script (test.py): Upper=+Γ, Lower=-Γ
   - **Both are correct** depending on convention, but must be consistent

3. **Position:**
   - Dynamic: Shed from separation point on cylinder surface
   - Static: Placed at geometric grid positions

---

## 7. SELF-INDUCTION PREVENTION

### Problem:
A vortex would induce **infinite velocity** on itself if treated as point vortex:
```
v_self = Γ/(2πr) → ∞ as r → 0
```

### Solutions:

| Script | Method | Cutoff Value | Implementation |
|--------|--------|-------------|----------------|
| **von_karman** | Exclude index | ε = 0.01D | Pass `exclude_self_idx` to velocity |
| **simple_vortex** | Distance check | 0.02D | `if dist < 0.02*D: continue` |
| **test.py** | Distance check | 1e-3 m | `if r < 1e-3: continue` |
| **test01.py** | Distance check | 1e-6 m | `if r < 1e-6: return 0.0, 0.0` |

### Code Examples:

**von_karman (index exclusion):**
```python
def compute_velocity_at_point(x, y, vortices, t, params, exclude_idx=-1):
    u_total = params['U_inf']
    v_total = 0.0
    
    for i, vortex in enumerate(vortices):
        if i == exclude_idx:  # Skip self
            continue
        
        # Compute contribution
        u, v = gaussian_velocity_contribution_dynamic(...)
        u_total += u
        v_total += v
    
    return u_total, v_total
```

**simple_vortex (distance check):**
```python
for xv, yv, gamma, t_birth in vortices:
    dist = np.sqrt((x - xv)**2 + (y - yv)**2)
    if dist < 0.02 * D:  # Self-induction cutoff
        continue
    
    # Compute contribution
    ...
```

**test01 (early return):**
```python
def induced_velocity(dx, dy, r, gamma, sigma):
    if r < 1e-6:  # Self cutoff
        return 0.0, 0.0
    
    # Compute contribution
    ...
```

### Physics:
- **Reality**: Vortex has finite core (σ), no true singularity
- **Lamb-Oseen model**: Regularizes singularity: V_θ → 0 as r → 0
- **Numerical**: Still need cutoff to avoid numerical issues

---

## 8. WAKE SPREADING (von_karman only)

### Mathematical Model:

**Transverse Spacing Growth:**
```
h(x) = h₀ + α·x

where:
  h₀ = initial spacing (von Kármán value: 0.281·D)
  α = spreading coefficient (typically 0.01-0.10)
  x = downstream distance
```

**Core Radius Growth (spatial):**
```
σ(x) = √[σ₀² + 4ν·(x/U∞)]
```

### Implementation (von_karman only):

```python
def compute_transverse_spacing(x, params):
    """
    Wake spreading model: vortex rows diverge downstream
    """
    h_0 = params['h']  # Initial spacing
    alpha = params['spreading_coef']  # Spreading rate
    return h_0 + alpha * x

def compute_vortex_core_radius(x, params):
    """
    Spatial core growth (for algebraic positioning)
    """
    sigma_0 = params['sigma_0']
    nu = params['U_inf'] * params['D'] / params['Re']
    return np.sqrt(sigma_0**2 + 4*nu*x/params['U_inf'])

# When positioning vortices:
y_pos = np.array([compute_transverse_spacing(x, params)/2 for x in x_pos])
sigma_pos = np.array([compute_vortex_core_radius(x, params) for x in x_pos])
```

### Visual Effect:
```
Without spreading (h = const):
  +  +  +  +  +  +  +  +     Upper row (parallel)
  -  -  -  -  -  -  -  -     Lower row (parallel)

With spreading [h(x) = h₀ + α·x]:
  +   +    +     +      +    Upper row (diverging)
  -   -    -     -      -    Lower row (diverging)
```

### Other Scripts:
- **simple_vortex**: Uses constant h, temporal σ(t)
- **test.py**: Uses constant h, spatial σ(x) but no α·x term
- **test01**: Uses constant h, temporal σ(t)

**Only von_karman** includes wake spreading model!

---

## 9. VELOCITY FIELD COMPUTATION

### Total Velocity at Point (x,y):

```
V_total = V_∞ + Σ V_vortex,i + Σ V_image,j

where:
  V_∞ = freestream (U∞, 0)
  V_vortex,i = contribution from real vortex i
  V_image,j = contribution from image vortex j (if BC enabled)
```

### Computational Complexity:

| Script | Real Vortices | Image Vortices | Total Cost |
|--------|--------------|----------------|------------|
| **von_karman** | N | 0 or N | O(N²) or O(2N²) |
| **simple_vortex** | N | 0 | O(N²) |
| **test.py** | N | 0 | O(N²) |
| **test01** | N | 4N (4 cylinders) | O(5N²) |

For test01 with 100 vortices:
- Real: 100 vortices
- Images: 400 images (100 vortices × 4 cylinders)
- Total: 500 contributions per evaluation point

---

## 10. NUMERICAL PARAMETERS

### Typical Values:

| Parameter | von_karman | simple_vortex | test.py | test01 |
|-----------|-----------|---------------|---------|--------|
| **D (diameter)** | 1.0 m | 1.0 m | 5.0 m | 4.0 m |
| **U∞ (velocity)** | 1.0 m/s | 1.5 m/s | 1.0 m/s | 1.0 m/s |
| **Re (Reynolds)** | 60 | 150 | 100 | 100 |
| **St (Strouhal)** | 0.2 | 0.2 | 0.2 | 0.2 |
| **σ₀ (core)** | 0.1D | 0.1D | 0.1D | 0.1D |
| **dt (timestep)** | 0.02 | 0.015 | N/A | 0.01 |
| **θ_sep (angle)** | 1.4 rad | 70° | N/A | 1.4 rad |

### Timestep Requirements:

**CFL Condition:**
```
dt < σ₀/U_max

For σ₀ = 0.1D, U_max ≈ 2U∞:
  dt < 0.05·(D/U∞)

With D=1, U∞=1:
  dt < 0.05
```

### Chosen Values:
- **von_karman**: dt = 0.02 ✓
- **simple_vortex**: dt = 0.015 ✓
- **test01**: dt = 0.01 ✓ (most conservative)

---

## SUMMARY

### Core Algorithm Differences:

1. **Time Integration:**
   - RK4 (test01): O(dt⁴) accuracy
   - Euler (von_karman, simple_vortex): O(dt) accuracy
   - None (test.py): Static

2. **Boundary Conditions:**
   - Image vortices (von_karman, test01): Enforced
   - None (simple_vortex, test.py): Violated

3. **Core Growth:**
   - Temporal σ(t) (von_karman, simple_vortex, test01): Correct
   - Spatial σ(x) (test.py): Approximate

4. **Wake Model:**
   - Spreading h(x) (von_karman only): Advanced
   - Constant h (others): Standard

### Most Accurate: **test01.py**
- RK4 integration
- Image vortices
- Multi-cylinder capability
- Proper temporal evolution

### Best for Education: **simple_vortex_street.py**
- Clear dynamic evolution
- Fast computation
- Easy to understand

### Most Flexible: **von_karman_simulator.py**
- Optional image vortices
- Wake spreading model
- Hybrid approach

### Least Accurate: **test.py**
- No dynamics
- No boundary conditions
- Kinematic approximation only
