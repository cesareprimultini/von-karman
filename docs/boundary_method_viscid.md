  Key Findings

  1. Your Image Vortex Implementation is Mathematically 
  Correct... for Inviscid Flow

  Looking at rk4-vortex-fast.py:78-90, your calculation
  is correct:
  - Position: x_img = x_cyl + (a²·x_rel) / r² ✓
  - Strength: gamma_img = -gamma ✓
  - This matches the classical method of images from
  potential flow theory

  2. But Image Vortex Method is for Inviscid Potential 
  Flow Only

  From the literature
  (https://open.oregonstate.education/intermediate-fluid-
  mechanics/chapter/the-panel-method-an-introduction/):

  "The method, as we present it here, is based on 
  inviscid flow analysis... cannot directly account for 
  viscous friction forces... assumes potential flow 
  satisfying Laplace's equation for the velocity 
  potential, which excludes viscous effects by 
  definition."

  The http://web.mit.edu/fluids-modules/www/potential_flo
  ws/LecturesHTML/lec1011/node37.html and https://farside
  .ph.utexas.edu/teaching/336L/Fluidhtml/node86.html
  confirm this applies strictly to potential flow.

  3. Modern Viscous Vortex Methods DON'T Use Image 
  Vortices

  The authoritative work by Koumoutsakos & Leonard (1995)
   on https://www.cambridge.org/core/journals/journal-of-
  fluid-mechanics/article/abs/highresolution-simulations-
  of-the-flow-around-an-impulsively-started-cylinder-usin
  g-vortex-methods/615FA46CB8A5BACF14FB4100A3CBC598 and
  their 1994 paper https://www.semanticscholar.org/paper/
  Boundary-Conditions-for-Viscous-Vortex-Methods-Koumouts
  akos-Leonard/adf58a03207c48783f7f6bb6710fbeda82129399
  states:

  "The vorticity creation process at the boundary, due to
   the no-slip condition, is expressed in terms of a 
  vorticity flux. The no-slip condition is not enforced 
  by the generation of new vortices at the boundary but 
  instead by modifying the strength of the vortices in 
  the vicinity of the boundary."

  They explicitly avoid using image vortices for viscous
  flow!

  4. Why You're Seeing Far Upstream Influence

  Several compounding factors:

  a) Cumulative Effect of Many Vortices:
  - At t=300s, you have potentially hundreds of vortices
  - Each has an image vortex (so double the count)
  - Even though individual vortex influence decays as
  1/r, with N vortices + N images = 2N sources,
  cumulative effect is significant
  - 4 cylinders multiplies this by 4

  b) Non-Physical Boundary Condition Enforcement:
  - Image vortices in potential flow create zero normal
  velocity at cylinder surface
  - But they don't enforce no-slip (zero tangential
  velocity)
  - Your viscous vortices are shedding and diffusing, but
   the boundary condition treatment is inviscid
  - This creates an inconsistent flow field that
  manifests as non-physical upstream influence

  c) The Inconsistency Propagates Upstream:
  From literature on https://open.oregonstate.education/i
  ntermediate-fluid-mechanics/chapter/potential-flows/:
  "Potential flow is unable to satisfy required boundary 
  conditions, especially near solid boundaries, which 
  makes it invalid in representing the required flow 
  field."

  Code-Specific Issues

  Lines 78-107 (during advection):

  # You compute image vortices for ALL vortices at once
  for cyl in cylinders:
      # For EACH cylinder, you create images of ALL
  vortices
      # This is correct for inviscid flow but problematic
   for viscous

  Lines 229-254 (velocity field computation):

  Same issue - you're creating image vortices assuming
  inviscid conditions.

  What Should You Do?

  Option 1: Acknowledge the Limitation (Simple)

  Accept that this is a simplified model mixing inviscid
  boundary treatment with viscous vortex evolution.
  Document this limitation but recognize it's
  computationally cheap and gives qualitatively correct
  wake patterns.

  Option 2: Remove Image Vortices (Moderate)

  - Remove all image vortex code (lines 77-107, 229-254)
  - Accept that vortices will penetrate the cylinders
  slightly
  - This is actually more consistent with your viscous
  approach
  - Use smaller time steps and rely on the shedding
  mechanism to create the wake

  Option 3: Implement Proper Viscous Boundary Conditions 
  (Advanced)

  Implement the vorticity flux method from Koumoutsakos &
   Leonard:
  - Calculate slip velocity at cylinder surface
  - Create new vorticity at the boundary based on viscous
   flux: ω_flux = -2ν(∂u_slip/∂n)
  - Modify nearby vortex strengths rather than using
  image vortices
  - Much more complex, but physically correct

  Option 4: Panel Method for Boundary (Hybrid)

  Replace scattered vortices near cylinders with panel 
  methods on the cylinder surface, combine with free
  vortices in wake. This is common in aerodynamics codes.

  My Recommendation

  For your application, I'd suggest Option 2 (remove
  image vortices) or Option 1 (document the limitation).
  Here's why:

  1. Your current approach is a hybrid that's neither
  fully inviscid nor fully viscous
  2. The far-upstream influence suggests the image
  vortices are creating non-physical flow patterns
  3. Simply letting vortices shed naturally (without
  image vortices) would be more consistent with your
  Lamb-Oseen viscous core model
  4. Many successful Von Kármán simulations use vorticity
   creation at the boundary rather than image vortices

  Additional Observations

  Your other parameters look reasonable:
  - sigma_0 = 0.1*D ✓
  - Gamma_mag = 2π·St·U_inf·D ✓
  - theta_sep = 80° ✓
  - Viscous core growth: sigma² = sigma_0² + 4νt ✓ (line
  149)

  Sources

  - https://open.oregonstate.education/intermediate-fluid
  -mechanics/chapter/the-panel-method-an-introduction/
  - http://web.mit.edu/fluids-modules/www/potential_flows
  /LecturesHTML/lec1011/node37.html
  - https://farside.ph.utexas.edu/teaching/336L/Fluidhtml
  /node86.html
  - https://www.cambridge.org/core/journals/journal-of-fl
  uid-mechanics/article/abs/highresolution-simulations-of
  -the-flow-around-an-impulsively-started-cylinder-using-
  vortex-methods/615FA46CB8A5BACF14FB4100A3CBC598
  - https://www.semanticscholar.org/paper/Boundary-Condit
  ions-for-Viscous-Vortex-Methods-Koumoutsakos-Leonard/ad
  f58a03207c48783f7f6bb6710fbeda82129399
  - https://www.mdpi.com/2311-5521/6/2/68
  - https://github.com/aGfrost/2D_Vortex_Methods
  - https://www.cambridge.org/core/journals/journal-of-fl
  uid-mechanics/article/vorticity-dynamics-at-partialslip
  -boundaries/84654172D32336F063F3191CDCCFE084