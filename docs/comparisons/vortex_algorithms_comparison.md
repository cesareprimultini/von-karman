# Von Kármán Vortex Street: Comprehensive Algorithm Comparison

## EXECUTIVE SUMMARY

Four different approaches to simulating von Kármán vortex streets, ranging from static algebraic models to fully dynamic time-stepping simulations:

**Algorithm Classification:**
1. **von_karman_simulator.py**: HYBRID (algebraic + dynamic) - Single cylinder
2. **simple_vortex_street.py**: DYNAMIC EULER - Single cylinder, no boundary conditions
3. **test.py**: STATIC ALGEBRAIC - 4 cylinders, no time evolution
4. **test01.py**: DYNAMIC RK4 - 4 cylinders, full boundary conditions

**Physical Realism Ranking (Low → High):**
```
test.py < von_karman_simulator.py < simple_vortex_street.py < test01.py
  (static)    (semi-dynamic)        (dynamic Euler)        (dynamic RK4 + BC)
```

---

## I. GENERAL OVERVIEW

### **von_karman_simulator.py** - HYBRID APPROACH
- **Type**: Algebraic positioning + Dynamic tracking option
- **Time Integration**: Euler method (when dynamic mode enabled)
- **Vortex Model**: Lamb-Oseen with viscous spreading
- **Boundary Conditions**: Image vortices (optional)
- **Key Feature**: Pre-positioned vortices at algebraic locations with optional dynamic evolution
- **Cylinders**: Single
- **Physical Realism**: Medium (primarily algebraic, can do dynamic)

### **simple_vortex_street.py** - PURE DYNAMIC SIMULATION  
- **Type**: Time-dependent vortex shedding
- **Time Integration**: Euler method
- **Vortex Model**: Lamb-Oseen with time-dependent core growth
- **Boundary Conditions**: NONE (vortices shed from surface, no image vortices)
- **Key Feature**: Vortices shed periodically from cylinder, advected by velocity field
- **Cylinders**: Single
- **Physical Realism**: Medium-High (dynamic but no proper BC)

### **test.py** - STATIC ALGEBRAIC
- **Type**: Steady-state visualization only
- **Time Integration**: NONE (no time evolution)
- **Vortex Model**: Lamb-Oseen with age-based core growth
- **Boundary Conditions**: NONE
- **Key Feature**: Fixed vortex positions assuming convection at U∞
- **Cylinders**: 4 (grid pattern)
- **Physical Realism**: Low (purely kinematic, no dynamics)

### **test01.py** - FULL DYNAMIC SIMULATION
- **Type**: Time-dependent vortex shedding with proper boundary conditions
- **Time Integration**: RK4 (4th order Runge-Kutta)
- **Vortex Model**: Lamb-Oseen with time-dependent core growth
- **Boundary Conditions**: Image vortices (enforces no-penetration)
- **Key Feature**: Most complete physical model with accurate time-stepping
- **Cylinders**: 4 (grid pattern)
- **Physical Realism**: HIGH (accurate integration + proper BC)

---

## II. DETAILED FORMULA-BY-FORMULA COMPARISON

### A. VORTEX CORE MODEL (σ - Core Radius)

All scripts use the **Lamb-Oseen vortex model**, but compute core radius differently:

| Script | Formula | Implementation | Notes |
|--------|---------|----------------|-------|
| **von_karman** | σ²(t) = σ₀² + 4νt | `sigma_sq = params['sigma_0']**2 + 4 * nu * age` | Time-dependent, age = t - t_birth |
| **simple_vortex** | σ²(t) = σ₀² + 4νt | `sigma = np.sqrt(sigma_0**2 + 4 * nu * age)` | Identical implementation |
| **test.py** | σ²(x) = σ₀² + 4ν(x/U∞) | `age_up = (x_up - x_cyl) / U_inf` | **Spatial age** (assumes U∞ convection) |
| **test01.py** | σ²(t) = σ₀² + 4νt | `v['sigma'] = np.sqrt(sigma_0**2 + 4 * nu * age)` | Time-dependent, updated each step |

**Key Difference**: test.py uses **spatial** age (x/U∞) while others use **temporal** age (t - t_birth)
- test.py: `age = distance/velocity` (kinematic assumption)
- Others: `age = current_time - birth_time` (true dynamic tracking)

**Physics**: 
- σ₀ = initial core radius (typically 0.1D)
- ν = kinematic viscosity = U∞·D/Re
- Core growth captures viscous diffusion

---

### B. LAMB-OSEEN VELOCITY FIELD

All scripts compute velocity using the Lamb-Oseen formula:

**Mathematical Formula:**
```
V_θ(r) = (Γ/2πr) × [1 - exp(-r²/σ²)]

u = -V_θ × (dy/r)  [x-component]
v = +V_θ × (dx/r)  [y-component]
```

Where:
- r² = dx² + dy² (distance from vortex center)
- Γ = circulation (signed)
- σ = core radius

#### Implementation Comparison:

**von_karman_simulator.py:**
```python
def gaussian_velocity_contribution(X, Y, x_v, y_v, sigma, gamma):
    dx = X - x_v
    dy = Y - y_v
    r_sq = dx**2 + dy**2
    r = np.sqrt(r_sq + 1e-10)
    
    V_theta = (gamma / (2 * np.pi * r)) * (1 - np.exp(-r_sq / (sigma**2 + 1e-10)))
    
    U = -V_theta * dy / r
    V = V_theta * dx / r
    return U, V
```

**simple_vortex_street.py:**
```python
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
```

**test.py:**
```python
sigma = vortex['sigma']
exp_factor = np.exp(-r**2 / sigma**2)
V_theta = (vortex['gamma'] / (2 * np.pi * r)) * (1 - exp_factor)

u += -V_theta * dy / r
v += +V_theta * dx / r
```

**test01.py:**
```python
def induced_velocity(dx, dy, r, gamma, sigma):
    if r < 1e-6:
        return 0.0, 0.0
    exp_factor = np.exp(-r**2 / sigma**2)
    V_theta = (gamma / (2 * np.pi * r)) * (1 - exp_factor)
    u_ind = -V_theta * dy / r
    v_ind = V_theta * dx / r
    return u_ind, v_ind
```

**Differences:**
- **Numerical stability**: All add small epsilon (1e-10 or 1e-6) to prevent division by zero
- **Vectorization**: von_karman uses full NumPy arrays, others compute point-by-point or vectorize later
- **Result**: IDENTICAL mathematical formulation

---

### C. CIRCULATION MAGNITUDE (Γ)

Different formulas for vortex strength:

| Script | Formula | Value | Physical Basis |
|--------|---------|-------|----------------|
| **von_karman** | Γ = C_γ × St × U∞ × D | C_γ = 2π (tunable) | From Strouhal scaling |
| **simple_vortex** | Γ = 2.5 × U∞ × D | Fixed constant | Empirical |
| **test.py** | Γ = 2π × St × U∞ × D | St = 0.2 | Standard formula |
| **test01.py** | Γ = 2π × St × U∞ × D | St = 0.2 | Standard formula |

**Standard Physics:**
- Strouhal number: St = f·D/U∞ ≈ 0.2 for Re ~ 100
- Shedding frequency: f = St·U∞/D
- Circulation: Γ ~ St·U∞·D (dimensionally consistent)

**Differences:**
- simple_vortex uses **empirical constant 2.5** (not physics-based)
- Others use **2π×St** which comes from Γ = ∮v·dl for a vortex

---

### D. VORTEX ADVECTION (Time Integration)

This is where the **biggest differences** occur:

#### **1. von_karman_simulator.py - EULER METHOD**
```python
def advect_vortices_euler(vortices, params, t, dt):
    new_vortices = []
    for vortex in vortices:
        u, v = compute_velocity_at_point(vortex.x, vortex.y, vortices, t, params)
        
        x_new = vortex.x + u * dt
        y_new = vortex.y + v * dt
        
        new_vortices.append(Vortex(x_new, y_new, vortex.gamma, vortex.t_birth))
    return new_vortices
```

**Order**: First-order (O(dt))
**Accuracy**: Low for large dt
**Stability**: Requires small dt (typically < 0.02)

#### **2. simple_vortex_street.py - EULER METHOD**
```python
for xv, yv, gamma, t_birth in vortices:
    u, v = get_velocity_at_point(xv, yv, vortices, t)
    
    xv_new = xv + u * dt
    yv_new = yv + v * dt
    
    if xv_new < 15*D and abs(yv_new) < 4*D:
        new_vortices.append([xv_new, yv_new, gamma, t_birth])
```

**Order**: First-order (O(dt))
**Same as**: von_karman_simulator.py
**Key Addition**: Domain boundary checks to remove vortices

#### **3. test.py - NO ADVECTION**
```python
# Static positions, no time evolution
x_up = x_start + k * l
y_up = y_cyl + h/2
```

**Order**: N/A (no dynamics)
**Assumption**: Vortices convect at U∞ in straight lines

#### **4. test01.py - RK4 METHOD (Most Accurate)**
```python
def rk4_advect(vortex_idx, vortices, dt):
    def vel_func(pos):
        return velocity_at_point(pos[0], pos[1], vortices, 
                                include_images=True, 
                                exclude_self_idx=vortex_idx)
    
    x0, y0 = vortices[vortex_idx]['x'], vortices[vortex_idx]['y']
    pos0 = np.array([x0, y0])
    
    k1 = np.array(vel_func(pos0)) * dt
    k2 = np.array(vel_func(pos0 + 0.5 * k1)) * dt
    k3 = np.array(vel_func(pos0 + 0.5 * k2)) * dt
    k4 = np.array(vel_func(pos0 + k3)) * dt
    
    pos_new = pos0 + (k1 + 2*k2 + 2*k3 + k4) / 6
    return pos_new[0], pos_new[1]
```

**Order**: Fourth-order (O(dt⁴))
**Accuracy**: HIGH - can use larger dt than Euler
**Key Feature**: Four velocity evaluations per step

**RK4 vs Euler Comparison:**
```
Euler:   x_{n+1} = x_n + v(x_n)·dt
RK4:     x_{n+1} = x_n + (k1 + 2k2 + 2k3 + k4)/6

Error per step:
  Euler: O(dt²)
  RK4:   O(dt⁵)

Total error after T:
  Euler: O(dt)
  RK4:   O(dt⁴)
```

**Impact**: RK4 can be **16× more accurate** for same dt, or use **2× larger dt** for same accuracy

---

### E. IMAGE VORTICES (Boundary Conditions)

Critical for enforcing no-penetration condition on cylinder surface.

#### **Mathematical Theory:**
For a cylinder of radius `a` centered at origin, a vortex at (x,y) creates an image at:
```
z_img = a²/z̄   (complex notation)

In Cartesian:
x_img = a²·x/(x² + y²)
y_img = a²·y/(x² + y²)
Γ_img = -Γ
```

| Script | Image Vortices | Implementation | Accuracy |
|--------|---------------|----------------|----------|
| **von_karman** | **YES (optional)** | Full implementation with validity check | HIGH |
| **simple_vortex** | **NO** | Not implemented | POOR (no BC) |
| **test.py** | **NO** | Not implemented | POOR (no BC) |
| **test01.py** | **YES** | Full implementation with validity check | HIGH |

#### **von_karman_simulator.py:**
```python
def get_image_vortex_position(vortex, params):
    x, y = vortex.x, vortex.y
    a = params['D'] / 2
    r_sq = x**2 + y**2
    
    if r_sq <= a**2:  # Vortex inside cylinder (invalid)
        return None
    
    x_img = (a**2 * x) / r_sq
    y_img = (a**2 * y) / r_sq
    return (x_img, y_img)

# Then in velocity calculation:
if params.get('use_image_vortices', False):
    image_pos = get_image_vortex_position(vortex, params)
    if image_pos is not None:
        U_img, V_img = gaussian_velocity_contribution(
            X, Y, x_img, y_img, sigma, -vortex.gamma  # Opposite sign!
        )
        U += U_img
        V += V_img
```

#### **test01.py:**
```python
def get_image_vortex(vortex, cyl):
    x_rel = vortex['x'] - cyl['x']
    y_rel = vortex['y'] - cyl['y']
    r_sq = x_rel**2 + y_rel**2
    
    if r_sq <= a**2 + 1e-6:
        return None
    
    x_img_rel = (a**2 * x_rel) / r_sq
    y_img_rel = (a**2 * y_rel) / r_sq
    x_img = cyl['x'] + x_img_rel
    y_img = cyl['y'] + y_img_rel
    gamma_img = -vortex['gamma']  # Opposite circulation
    sigma_img = vortex['sigma']
    
    return {'x': x_img, 'y': y_img, 'gamma': gamma_img, 'sigma': sigma_img}
```

**Key Differences:**
1. test01.py handles **multiple cylinders** (4 cylinders, 4 images per vortex)
2. test01.py uses **cylinder-relative coordinates** (x_rel, y_rel)
3. Both correctly implement **opposite circulation** (Γ_img = -Γ)

**Physical Impact:**
- **Without images**: Vortices can pass through cylinder (non-physical)
- **With images**: No-penetration BC enforced, vortices deflect around cylinder

---

### F. VORTEX SHEDDING

How new vortices are introduced:

| Script | Shedding Method | Location | Sign Convention |
|--------|----------------|----------|-----------------|
| **von_karman** | Periodic at t = n·T_shed | θ = π ± θ_sep | Upper: -Γ (CW), Lower: +Γ (CCW) |
| **simple_vortex** | Periodic at t = n·T_shed | θ = π ± θ_sep | Upper: -Γ (CW), Lower: +Γ (CCW) |
| **test.py** | N/A (pre-positioned) | Algebraic grid | Upper: +Γ, Lower: -Γ |
| **test01.py** | Periodic at t = n·T_shed | θ = π ± θ_sep | Upper: -Γ (CW), Lower: +Γ (CCW) |

#### **Shedding Period:**
```
T_shed = D / (St × U∞)
```

Where St ≈ 0.2 for Re ~ 100-200.

#### **Separation Angle:**
- **Typical value**: θ_sep ≈ 70-80° from rear stagnation point
- Physically: Flow separates at this angle on cylinder surface

#### **von_karman & simple_vortex & test01:**
```python
T_shed = D / (St * U_inf)
expected_count = int(t / T_shed)

while shed_count < expected_count:
    is_upper = (shed_count % 2 == 0)  # Alternate rows
    
    theta = np.pi - theta_sep if is_upper else np.pi + theta_sep
    x_shed = (D/2 + sigma_0) * np.cos(theta)
    y_shed = (D/2 + sigma_0) * np.sin(theta)
    
    gamma = -Gamma_shed if is_upper else Gamma_shed
    
    vortices.append([x_shed, y_shed, gamma, t])
    shed_count += 1
```

**Key Points:**
- Shed **just outside cylinder surface** (radius = D/2 + σ₀)
- **Alternating pattern**: upper → lower → upper → lower
- **Opposite circulation**: ensures alternating street pattern

#### **test.py (algebraic):**
```python
for k in range(n_pairs):
    # Upper row: positive circulation
    x_up = x_start + k * l
    y_up = y_cyl + h/2
    vortices.append({'x': x_up, 'y': y_up, 'gamma': +Gamma, ...})
    
    # Lower row: negative circulation  
    x_low = x_start + a + k * l
    y_low = y_cyl - h/2
    vortices.append({'x': x_low, 'y': y_low, 'gamma': -Gamma, ...})
```

**Difference**: Pre-places ALL vortices at **fixed geometric positions**
- No time evolution
- Assumes steady-state convection at U∞
- Uses von Kármán stability criterion: h/l ≈ 0.281

---

### G. SELF-INDUCTION CUTOFF

Prevents vortices from inducing infinite velocity on themselves:

| Script | Cutoff Radius | Implementation |
|--------|---------------|----------------|
| **von_karman** | ε = 0.01D | `exclude_self_idx` in RK4 |
| **simple_vortex** | 0.02D | `if dist < 0.02 * D: continue` |
| **test.py** | 1e-3 m | `if r < 1e-3: continue` |
| **test01.py** | 1e-6 m | `if r < 1e-6: return 0.0, 0.0` |

**Methods:**
1. **Skip in velocity loop**: When computing total velocity, skip self-contribution
2. **Exclude from RK4**: Pass vortex index to exclude from velocity calculation
3. **Cutoff radius**: If distance < ε, set velocity to zero

**Physics**: In reality, vortex has finite core - singularity is regularized by σ

---

### H. WAKE SPREADING (von_karman only)

**Unique to von_karman_simulator.py:**

```python
def compute_transverse_spacing(x, params):
    """Wake spreading: h(x) = h₀ + α·x"""
    h_0 = params['h']
    alpha = params['spreading_coef']
    return h_0 + alpha * x

def compute_vortex_core_radius(x, params):
    """Core growth: σ(x) = √(σ₀² + 4ν·x/U∞)"""
    sigma_0 = params['sigma_0']
    nu = params['U_inf'] * params['D'] / params['Re']
    return np.sqrt(sigma_0**2 + 4 * nu * x / params['U_inf'])
```

**Physical Meaning:**
- Vortex rows **diverge** downstream (not parallel)
- Core radius **grows** with downstream distance
- Models **viscous spreading** of wake

**Not in other scripts** - they use **constant h** and **time-based σ(t)**

---

## III. SUMMARY TABLE: ALGORITHM COMPARISON

| Feature | von_karman | simple_vortex | test.py | test01.py |
|---------|------------|---------------|---------|-----------|
| **Type** | Hybrid | Dynamic | Static | Dynamic |
| **Time Integration** | Euler (optional) | Euler | None | RK4 |
| **Integration Order** | O(dt) | O(dt) | N/A | O(dt⁴) |
| **Vortex Model** | Lamb-Oseen | Lamb-Oseen | Lamb-Oseen | Lamb-Oseen |
| **Core Growth** | σ(t) or σ(x) | σ(t) | σ(x) | σ(t) |
| **Wake Spreading** | YES (h(x)) | NO | NO | NO |
| **Image Vortices** | Optional | NO | NO | YES |
| **Cylinders** | 1 | 1 | 4 | 4 |
| **Circulation** | 2π·St·U∞·D | 2.5·U∞·D | 2π·St·U∞·D | 2π·St·U∞·D |
| **Shedding** | Dynamic | Dynamic | Algebraic | Dynamic |
| **Boundary Conditions** | Image vortices | None | None | Image vortices |
| **Self-Induction** | Exclude index | Distance cutoff | Distance cutoff | Distance cutoff |
| **Computational Cost** | Medium | Low | Very Low | HIGH |
| **Physical Realism** | Medium | Medium | Low | **HIGHEST** |

---

## IV. ACCURACY AND PHYSICAL REALISM ANALYSIS

### **Ranking by Physical Accuracy:**

#### **1. test01.py - HIGHEST**
**Strengths:**
- ✅ RK4 integration (O(dt⁴) accuracy)
- ✅ Image vortices (proper boundary conditions)
- ✅ Time-dependent core growth
- ✅ Self-exclusion in velocity calculation
- ✅ Multiple cylinder interactions

**Weaknesses:**
- ⚠️ High computational cost (4× evaluations per step × N_vortices × N_cylinders)
- ⚠️ No viscous spreading (h constant)

**Best for**: Accurate multi-cylinder wake interaction studies

---

#### **2. simple_vortex_street.py - MEDIUM-HIGH**
**Strengths:**
- ✅ True dynamic shedding
- ✅ Time-dependent core growth
- ✅ Simple and efficient

**Weaknesses:**
- ❌ No image vortices (boundary conditions violated)
- ❌ Euler integration (O(dt) accuracy)
- ❌ Empirical circulation formula (not physics-based)

**Best for**: Quick visualizations, educational demonstrations

---

#### **3. von_karman_simulator.py - MEDIUM**
**Strengths:**
- ✅ Optional image vortices
- ✅ Wake spreading model (viscous effects)
- ✅ Flexible (can run static or dynamic)

**Weaknesses:**
- ❌ Euler integration when dynamic
- ⚠️ Hybrid approach (mixing algebraic and dynamic)
- ⚠️ Complexity without clear advantage

**Best for**: Testing different models, exploring wake spreading

---

#### **4. test.py - LOWEST**
**Strengths:**
- ✅ Very fast (no time-stepping)
- ✅ Multi-cylinder capability

**Weaknesses:**
- ❌ No dynamics (purely kinematic)
- ❌ No boundary conditions
- ❌ Assumes steady-state immediately
- ❌ No vortex interactions (pre-positioned)

**Best for**: Quick static wake visualization only

---

## V. KEY ALGORITHMIC DIFFERENCES

### **A. Time-Stepping Method**

**Impact on Accuracy:**
```
Same dt = 0.01, same simulation time = 10s

Euler (simple_vortex, von_karman):
  - Error per step: O(dt²) = O(10⁻⁴)
  - Total error: O(dt) = O(10⁻²)
  - Requires ~1000 steps for decent accuracy

RK4 (test01):
  - Error per step: O(dt⁵) = O(10⁻¹⁰)
  - Total error: O(dt⁴) = O(10⁻⁸)
  - Can use dt = 0.04 for same accuracy (4× faster)
```

**Practical Impact:**
- RK4 captures **vortex pairing** and **complex interactions** accurately
- Euler may show **artificial drift** or **instability** at large dt
- For **long-time simulations**, RK4 essential

---

### **B. Boundary Conditions (Image Vortices)**

**Physical Impact:**

**Without Images (simple_vortex, test.py):**
- Vortices can **penetrate cylinder** (non-physical)
- No **deflection** around obstacles
- Wake develops as if **cylinder is permeable**

**With Images (test01, von_karman):**
- **No-penetration BC enforced**
- Vortices **deflect** realistically
- Velocity **exactly zero** normal to cylinder surface
- **Flow separation** captured correctly

**Example:**
```
At cylinder surface (r = D/2):
  Without images: v_normal ≠ 0 (flow through cylinder!)
  With images:    v_normal = 0 (physically correct)
```

---

### **C. Multi-Cylinder Interactions**

**test01.py**: Each vortex sees **4 image vortices** (one per cylinder)
```
Total velocity contributions per vortex:
  - 1× self (real vortex)
  - N-1× other real vortices
  - 4× image vortices (one per cylinder)
  
For N=100 vortices:
  - 99 real + 400 images = 499 interactions per vortex!
```

**Computational Cost:**
```
Simple (no images): O(N²)
With images:        O(N² + 4N) ≈ O(5N²) for 4 cylinders
```

---

## VI. RECOMMENDATIONS

### **For Educational/Visualization:**
→ **simple_vortex_street.py**
- Fast, clear, demonstrates key physics
- Good for understanding vortex shedding

### **For Research/Accurate Simulation:**
→ **test01.py**
- Most accurate integration
- Proper boundary conditions
- Handles multi-cylinder interactions

### **For Exploring Wake Spreading:**
→ **von_karman_simulator.py**
- Only one with h(x) spreading model
- Flexible testing environment

### **For Static Snapshots Only:**
→ **test.py**
- Very fast
- Limited physical meaning

---

## VII. POTENTIAL IMPROVEMENTS

### **For simple_vortex_street.py:**
1. Add image vortices → proper BC
2. Upgrade to RK4 integration
3. Use physics-based circulation (2π·St·U∞·D instead of 2.5·U∞·D)

### **For test01.py:**
1. Add wake spreading model h(x) = h₀ + α·x
2. Implement adaptive time-stepping (adjust dt based on local velocity)
3. Add vortex merging (when two vortices get very close)

### **For von_karman_simulator.py:**
1. Upgrade to RK4 when in dynamic mode
2. Clarify hybrid vs pure dynamic modes
3. Extend to multiple cylinders

---

## VIII. VALIDATION METRICS

All scripts should check these physics quantities:

### **A. Strouhal Number**
```
St = f·D/U∞ ≈ 0.2  (for Re ~ 100)

Measured from velocity probe time series
```

**Implemented in:**
- von_karman: ✅ (has probe time series)
- simple_vortex: ❌
- test.py: ❌ (no time evolution)
- test01: ❌

### **B. von Kármán Stability**
```
h/l ≈ 0.281  (stable configuration)
```

**Implemented in:**
- von_karman: ✅
- simple_vortex: Uses θ_sep instead
- test.py: ✅
- test01: Uses θ_sep instead

### **C. Vortex Convection Speed**
```
U_vortex ≈ U∞ + V₀  (where V₀ = orbital velocity)
```

Can measure from vortex trajectories in dynamic simulations.

---

## IX. CONCLUSION

**Most Physically Realistic:** test01.py (RK4 + images + proper BC)
**Most Educational:** simple_vortex_street.py (clear dynamics, fast)
**Most Flexible:** von_karman_simulator.py (hybrid approach, spreading model)
**Least Accurate:** test.py (static, no dynamics)

**Critical Formula Differences:**
1. **Time integration**: RK4 (O(dt⁴)) vs Euler (O(dt)) vs None
2. **Boundary conditions**: Image vortices (yes/no)
3. **Core growth**: σ(t) vs σ(x)
4. **Wake spreading**: h(x) vs constant h

For **accurate research simulations**, use **test01.py** or upgrade simple_vortex with RK4 and images.
