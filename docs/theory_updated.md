Here is a streamlined, mathematically precise explanation of how the velocity field is defined in this model. It bridges the gap between the concept and the specific equations needed to calculate velocity at any point $(x, y)$ and time $t$.

***

# The Time-Dependent Velocity Field

## The Core Concept: Superposition
The velocity at any point in the wake is calculated using the **Principle of Superposition**. The fluid flow is treated as a sum of two components:
1.  The constant background flow (Freestream).
2.  The induced rotation from every vortex currently in the simulation.

$$
\mathbf{u}(x,y,t) = \mathbf{U}_{\infty} + \sum_{i=1}^{N} \mathbf{u}_{i}(x,y,t)
$$

Where $\mathbf{u}_i$ is the velocity induced by the $i$-th vortex.

---

## 1. Vortex Generation: Periodic Shedding from the Cylinder

Before vortices can influence the flow field, they must first be created. In this model, vortices are **shed periodically** from the cylinder surface, alternating between the upper and lower separation points.

### Shedding Frequency
The shedding frequency is governed by the **Strouhal number** $St$, a dimensionless parameter that relates the vortex shedding frequency to the flow conditions:

$$
St = \frac{f \cdot D}{U_\infty}
$$

Where:
* **$f$:** Shedding frequency (vortices per second)
* **$D$:** Cylinder diameter
* **$U_\infty$:** Freestream velocity

For Reynolds numbers in the range $47 < Re < 190$, the Strouhal number is approximately $St \approx 0.2$, which yields a **shedding period**:

$$
T_{shed} = \frac{D}{St \cdot U_\infty}
$$

**Reference:** Roshko, A. (1954). *On the Development of Turbulent Wakes from Vortex Streets*. NACA Report 1191.

### Shedding Location and Circulation
Vortices are shed from the **separation points** on the cylinder surface, located at an angle $\theta_{sep}$ from the rear stagnation point (typically $\theta_{sep} \approx 80°$ or $1.4$ radians):

$$
\begin{aligned}
\text{Upper vortex position:} \quad & x = a\cos(\pi - \theta_{sep}), \quad y = a\sin(\pi - \theta_{sep}) \\
\text{Lower vortex position:} \quad & x = a\cos(\pi + \theta_{sep}), \quad y = a\sin(\pi + \theta_{sep})
\end{aligned}
$$

Where $a = D/2$ is the cylinder radius.

Each newly shed vortex carries a circulation magnitude:

$$
\Gamma_{shed} = C_\gamma \cdot St \cdot U_\infty \cdot D
$$

Where $C_\gamma$ is a tunable coefficient (typically $C_\gamma \approx 2\pi$ to $4\pi$). The vortices alternate in sign: **upper vortices are clockwise** (negative $\Gamma$) and **lower vortices are counter-clockwise** (positive $\Gamma$).

**Reference:** Williamson, C. H. K. (1996). *Vortex Dynamics in the Cylinder Wake*. Annual Review of Fluid Mechanics, 28, 477-539.

---

## 1.5. Gradual Circulation Buildup: The Formation Time

### The Problem with Instantaneous Vortex Creation

In a naive implementation, newly shed vortices appear with **full circulation strength** $\Gamma_{shed}$ instantaneously at $t = t_{birth}$. This creates a **discontinuity** in the velocity field: the vortex's induced velocity suddenly appears at all points in the domain, causing unphysical "jumps" in velocity measurements.

**Physical Reality:** Vortex shedding is a gradual process. Circulation builds up over time as the shear layer rolls up and detaches from the cylinder surface.

### The Formation Number Concept

Gharib, Rambod, and Shariff (1998) discovered experimentally that vortices reach their maximum circulation after a characteristic **formation time** defined by the dimensionless **formation number**:

$$
T^* = \frac{U_\infty \cdot t_{formation}}{D} \approx 3.6 - 4.5
$$

For most cylinder flows, $T^* \approx 4$ provides excellent agreement with experiments.

**Physical Interpretation:** A vortex requires approximately **4 diameters of travel time** (or equivalently, 4 convective time units) to fully develop its circulation.

**In dimensional form:**
$$
t_{formation} = \frac{4D}{U_\infty}
$$

For your simulation parameters ($D = 4.88$ m, $U_\infty = 1.5$ m/s), this gives $t_{formation} \approx 13$ seconds.

**Reference:** Gharib, M., Rambod, E., & Shariff, K. (1998). *A universal time scale for vortex ring formation*. Journal of Fluid Mechanics, 360, 121-140.

### Smooth Ramp Functions

To implement gradual circulation buildup, we apply a **smooth ramp function** $f(t)$ to the vortex circulation:

$$
\Gamma_{effective}(t) = \Gamma_{shed} \cdot f\left(\frac{t - t_{birth}}{t_{formation}}\right)
$$

Where $f(s)$ is a smooth function satisfying:
- $f(0) = 0$ (zero circulation at birth)
- $f(1) = 1$ (full circulation after formation time)
- $f'(s)$ is continuous (smooth velocity evolution)

#### Recommended: Smoothstep Function (Perlin, 1985)

The **smoothstep** (or "smootherstep") polynomial is widely used in computational fluid dynamics:

$$
f(s) = \begin{cases}
0 & \text{if } s \leq 0 \\
3s^2 - 2s^3 & \text{if } 0 < s < 1 \\
1 & \text{if } s \geq 1
\end{cases}
$$

**Properties:**
- $C^1$ continuous (smooth first derivative)
- Polynomial evaluation (computationally cheap)
- Exactly 0 before formation, exactly 1 after
- Symmetric acceleration/deceleration profile

**Implementation:**
```python
def smoothstep(x):
    """Perlin smoothstep function: 3x² - 2x³"""
    if x <= 0:
        return 0.0
    elif x >= 1:
        return 1.0
    else:
        return x * x * (3 - 2 * x)

# For each vortex:
age = t - t_birth
s = age / t_formation
gamma_effective = gamma * smoothstep(s)
```

**Reference:** Perlin, K. (1985). *An image synthesizer*. ACM SIGGRAPH Computer Graphics, 19(3), 287-296. (Smoothstep function widely adopted in computer graphics and CFD)

#### Alternative: Hyperbolic Tangent Ramp

For even smoother profiles ($C^\infty$ continuous), the hyperbolic tangent can be used:

$$
f(s) = \frac{1}{2}\left[1 + \tanh\left(4(s - 0.5)\right)\right]
$$

This provides infinite smoothness at the cost of slightly more computation (transcendental functions).

**Reference:** Standard practice in discrete vortex methods. See Cottet & Koumoutsakos (2000), *Vortex Methods*, Chapter 3.

#### Alternative: Raised Cosine (Hann Window)

From signal processing, provides smooth $C^1$ continuity:

$$
f(s) = \begin{cases}
0 & \text{if } s \leq 0 \\
\frac{1}{2}\left[1 - \cos(\pi s)\right] & \text{if } 0 < s < 1 \\
1 & \text{if } s \geq 1
\end{cases}
$$

### Numerical Implementation

In the velocity calculation routine, replace the raw circulation $\Gamma$ with the time-dependent effective circulation:

```python
for xv, yv, gamma, t_birth, cyl_idx in vortices:
    age = t - t_birth
    sigma = np.sqrt(sigma_0**2 + 4 * nu * age)

    # Apply smooth ramp to circulation
    s = age / t_formation
    gamma_effective = gamma * smoothstep(s)

    u, v = lamb_oseen_velocity(x, y, xv, yv, gamma_effective, sigma)
    u_total += u
    v_total += v
```

### Impact on Velocity Discontinuities

By applying the smooth ramp:
1. **Eliminates instantaneous appearance** of velocity contributions
2. **Continuous velocity field** in time (no jumps)
3. **More physically accurate** representation of the shear layer roll-up process
4. **Improves numerical stability** for velocity measurements near shedding regions

**Key Result:** Velocity measurements will show smooth oscillations corresponding to vortex shedding frequency, without artificial discontinuities when new vortices are created.

**Additional References:**
- **Rosenfeld, M., Rambod, E., & Gharib, M. (1998).** *Circulation and formation number of laminar vortex rings*. Journal of Fluid Mechanics, 376, 297-318.
- **Control and Suppression of Vortex Shedding** - MDPI Energies (2020). Demonstrates discrete vortex method implementations for cylinder flows. https://www.mdpi.com/1996-1073/13/17/4481
- **Vortex Shedding and Induced Forces** - Cambridge Core Review (2024). Modern perspective on computational vortex methods. https://www.cambridge.org/core/blog/2024/10/17/vortex-shedding-and-induced-forces-in-unsteady-flow/

---

## 2. The Building Block: The Lamb-Oseen Vortex
In this model, a vortex is not a singularity (a point of infinite spin). It is a **"Gaussian Blob"** of vorticity. This is the **Lamb-Oseen vortex**, which is an exact solution to the Navier-Stokes equations for a decaying vortex.

### The Magnitude Formula
The tangential (swirling) speed $V_\theta$ induced by a single vortex at a distance $r$ from its center is:

$$
V_\theta(r, t) = \frac{\Gamma}{2\pi r} \left[ 1 - \exp\left(-\frac{r^2}{\sigma(t)^2}\right) \right]
$$

* **$\Gamma$ (Gamma):** The circulation (strength) of the vortex. Positive for counter-clockwise, negative for clockwise.
* **$r$:** The distance from the probe point to the vortex center: $r = \sqrt{(x-x_i)^2 + (y-y_i)^2}$.
* **$\sigma(t)$:** The **Core Radius**. This is the critical viscous parameter.

> **Why this equation?**
> * **Far away ($r \gg \sigma$):** The exponential term vanishes. It simplifies to $\frac{\Gamma}{2\pi r}$, matching the classical inviscid point vortex.
> * **At the center ($r \to 0$):** The term in brackets approaches zero faster than $1/r$ blows up. The velocity goes to **zero**. This prevents the numerical "explosions" seen in inviscid models.

**Reference:** Lamb, H. (1932). *Hydrodynamics* (6th ed.). Cambridge University Press. (§§159-161: Vortex motion and diffusion).

---

## 3. Viscous Physics: Time-Dependent Core Growth
The "Viscous" part of the model appears in how the core size $\sigma$ changes. As time passes, viscosity causes the vorticity to diffuse (spread out), making the vortex larger and weaker at the center.

$$
\sigma^2(t) = \sigma_0^2 + 4\nu (t - t_{birth})
$$

* **$\sigma_0$:** The initial size of the vortex when shed from the cylinder (typically $\sigma_0 \approx 0.05D$ to $0.15D$).
* **$\nu$:** Kinematic viscosity ($\nu = U_\infty D / Re$).
* **$t - t_{birth}$:** The "age" of the vortex.

**Key Insight:** Older vortices downstream have larger $\sigma$ values. They are "fuzzier" and induce gentler velocities than new, sharp vortices near the cylinder.

**Reference:** Oseen, C. W. (1911). *Über Wirbelbewegung in einer reibenden Flüssigkeit*. Arkiv för Matematik, Astronomi och Fysik, 7(14).

---

## 4. Vortex Advection: Computing the Vortex Trajectories

Once shed, vortices do not remain stationary—they are **advected** (carried) by the local velocity field. The path of each vortex is computed by integrating the velocity at its position over time.

### The Advection Equation
Each vortex moves according to the ordinary differential equation (ODE):

$$
\frac{d\mathbf{x}_i}{dt} = \mathbf{u}(\mathbf{x}_i, t)
$$

Where $\mathbf{x}_i = (x_i, y_i)$ is the position of the $i$-th vortex, and $\mathbf{u}(\mathbf{x}_i, t)$ is the **total velocity field** at that location, computed using the superposition principle from Section 1 (excluding the vortex's own contribution to avoid singularities).

### Velocity Induced by Other Vortices
The velocity at vortex $i$ due to all other vortices $j \neq i$ is:

$$
\mathbf{u}(\mathbf{x}_i, t) = \mathbf{U}_\infty + \sum_{j \neq i} \frac{\Gamma_j}{2\pi r_{ij}} \left[1 - \exp\left(-\frac{r_{ij}^2}{\sigma_j^2}\right)\right] \hat{\mathbf{e}}_{\theta,ij}
$$

Where:
* $r_{ij} = |\mathbf{x}_i - \mathbf{x}_j|$ is the distance between vortices $i$ and $j$
* $\hat{\mathbf{e}}_{\theta,ij}$ is the tangential unit vector (perpendicular to the line connecting the vortices)

### Numerical Integration: Runge-Kutta 4th Order (RK4)
To solve the advection ODE, the model uses the **4th-order Runge-Kutta method** (RK4), which provides excellent accuracy with moderate computational cost:

$$
\begin{aligned}
\mathbf{k}_1 &= \mathbf{u}(t, \mathbf{x}_i) \\
\mathbf{k}_2 &= \mathbf{u}(t + \tfrac{\Delta t}{2}, \mathbf{x}_i + \tfrac{\Delta t}{2}\mathbf{k}_1) \\
\mathbf{k}_3 &= \mathbf{u}(t + \tfrac{\Delta t}{2}, \mathbf{x}_i + \tfrac{\Delta t}{2}\mathbf{k}_2) \\
\mathbf{k}_4 &= \mathbf{u}(t + \Delta t, \mathbf{x}_i + \Delta t\mathbf{k}_3) \\
\mathbf{x}_i(t + \Delta t) &= \mathbf{x}_i(t) + \frac{\Delta t}{6}(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4)
\end{aligned}
$$

**Stability Requirement:** The timestep must satisfy $\Delta t < \sigma_0 / U_{max}$ to prevent numerical instability. Typical values are $\Delta t \approx 0.02$ for $\sigma_0 = 0.1$ and $U_\infty = 1$.

**Reference:** Cottet, G. H., & Koumoutsakos, P. (2000). *Vortex Methods: Theory and Practice*. Cambridge University Press. (Chapter 2: Vortex advection schemes).

### Image Vortex Method for Cylinder Boundary Condition
To enforce the no-penetration boundary condition on the cylinder surface, an **image vortex** technique is employed. For each vortex at position $(x, y)$, an image vortex with **opposite circulation** $-\Gamma$ is placed at:

$$
x_{img} = \frac{a^2 x}{x^2 + y^2}, \quad y_{img} = -\frac{a^2 y}{x^2 + y^2}
$$

Where $a = D/2$ is the cylinder radius. The image vortices ensure that the net radial velocity at the cylinder surface is zero.

**Reference:** Milne-Thomson, L. M. (1968). *Theoretical Hydrodynamics* (5th ed.). Dover Publications. (Chapter 9: The method of images).

---

## 5. Cartesian Velocity Components (The Calculation Step)
To implement this in code (or finding the velocity at a specific point), you convert the tangential velocity $V_\theta$ into $u$ (x-velocity) and $v$ (y-velocity) components.

For a target point $(x, y)$ and a vortex located at $(x_i, y_i)$:

1.  **Calculate distance:** $r = \sqrt{(x - x_i)^2 + (y - y_i)^2}$
2.  **Calculate core size:** $\sigma = \sqrt{\sigma_0^2 + 4\nu \cdot \text{age}}$
3.  **Calculate Tangential Speed:** Use the Lamb-Oseen formula above to get $V_\theta$.
4.  **Project to Cartesian:**
    $$
    u_{induced} = -V_\theta \cdot \frac{(y - y_i)}{r}
    $$
    $$
    v_{induced} = +V_\theta \cdot \frac{(x - x_i)}{r}
    $$

*(Note: The negative sign in $u$ is because a counter-clockwise rotation at the top of a circle pushes fluid to the left)*.

**Reference:** Batchelor, G. K. (1967). *An Introduction to Fluid Dynamics*. Cambridge University Press. (§2.7: The velocity induced by a line vortex).

---

## 6. Key Assumptions and Limitations

### Assumptions
1.  **2D Flow:** The model assumes the cylinder is infinitely long and flow does not vary along the span ($z$-axis).
2.  **Linear Superposition:** We assume we can simply add the velocity fields of individual vortices together. This is strictly true for potential flow (inviscid) but is an approximation for viscous flow, acceptable because the vorticity is highly concentrated in the cores.
3.  **Passive Advection:** The vortices move at the local fluid velocity. They do not have their own "propulsion"; they simply drift in the flow created by all other vortices.
4.  **No Vortex Stretching:** The 2D assumption means vortex tubes cannot stretch or compress along their axes, which eliminates a key mechanism for energy cascade in 3D turbulence.

### Limitations
1.  **Near-Wake Inaccuracy:** In the "formation region" immediately behind the cylinder ($x < 1D$), the flow is a complex shear layer, not yet discrete vortices. This model approximates that region abruptly by shedding fully-formed vortices at the separation points.
2.  **No Vortex Merging:** In reality, two vortices might collide and merge. This algebraic model typically allows them to pass through or orbit each other without merging physics unless explicitly added.
3.  **Reynolds Number Range:** The model is most accurate for $47 < Re < 190$, where the wake is laminar and two-dimensional. Outside this range, three-dimensional effects or wake instabilities dominate.

---

## Scientific Sources

1.  **Lamb, H. (1932).** *Hydrodynamics*. Cambridge University Press. (The derivation of the viscous diffusion of a vortex line).
2.  **Oseen, C. W. (1911).** *Über Wirbelbewegung in einer reibenden Flüssigkeit*. Arkiv för Matematik, Astronomi och Fysik, 7(14). (Original Lamb-Oseen vortex solution).
3.  **Batchelor, G. K. (1967).** *An Introduction to Fluid Dynamics*. Cambridge University Press. (Standard reference for the superposition of line vortices and induced velocities).
4.  **Roshko, A. (1954).** *On the Development of Turbulent Wakes from Vortex Streets*. NACA Report 1191. (The physics of wake spreading, Reynolds number dependence, and Strouhal number measurements).
5.  **Williamson, C. H. K. (1996).** *Vortex Dynamics in the Cylinder Wake*. Annual Review of Fluid Mechanics, 28, 477-539. (Modern review of vortex shedding mechanisms and wake dynamics).
6.  **Cottet, G. H., & Koumoutsakos, P. (2000).** *Vortex Methods: Theory and Practice*. Cambridge University Press. (Modern implementation of Gaussian blobs for simulation and numerical advection schemes).
7.  **Milne-Thomson, L. M. (1968).** *Theoretical Hydrodynamics* (5th ed.). Dover Publications. (The method of images for boundary conditions in potential flow).
8.  **Gharib, M., Rambod, E., & Shariff, K. (1998).** *A universal time scale for vortex ring formation*. Journal of Fluid Mechanics, 360, 121-140. (Experimental discovery of the formation number for vortex circulation buildup).
9.  **Rosenfeld, M., Rambod, E., & Gharib, M. (1998).** *Circulation and formation number of laminar vortex rings*. Journal of Fluid Mechanics, 376, 297-318. (Detailed analysis of circulation growth during vortex formation).
10. **Perlin, K. (1985).** *An image synthesizer*. ACM SIGGRAPH Computer Graphics, 19(3), 287-296. (Introduction of the smoothstep function, widely adopted in CFD for smooth transitions).
11. **Control and Suppression of Vortex Shedding from a Slightly Rough Circular Cylinder by a Discrete Vortex Method.** MDPI Energies (2020). https://www.mdpi.com/1996-1073/13/17/4481 (Modern discrete vortex method implementation for cylinder flows).
12. **Vortex Shedding and Induced Forces in Unsteady Flow.** Cambridge Core Review (2024). https://www.cambridge.org/core/blog/2024/10/17/vortex-shedding-and-induced-forces-in-unsteady-flow/ (Contemporary perspective on computational vortex methods and validation).
