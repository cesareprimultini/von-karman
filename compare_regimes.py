"""
Create comparison plots showing laminar vs turbulent regime behavior
Demonstrates the effect of turbulence features on wake dynamics
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Key parameters
D = 4.88      # Cylinder diameter [m]
U_inf = 1.5   # Freestream velocity [m/s]

# Test cases
test_cases = [
    {'Re': 100, 'St': 0.167, 'regime': 'Laminar', 'color': 'blue'},
    {'Re': 6.7e6, 'St': 0.300, 'regime': 'Transcritical', 'color': 'red'}
]

def get_strouhal_number(Re):
    """Strouhal number as function of Reynolds number"""
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

def compute_eddy_viscosity_coeff(Re):
    """Eddy viscosity coefficient C_eddy as function of Re"""
    if Re < 1000:
        return 0.01
    elif Re < 1e5:
        return 0.02 + 0.01 * np.log10(Re / 1000)
    else:
        return 0.10

# Create comprehensive comparison figure
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

# ============================================================
# Plot 1: Strouhal Number vs Reynolds Number (Full Range)
# ============================================================
ax1 = fig.add_subplot(gs[0, :])

Re_range = np.logspace(2, 7, 200)
St_range = [get_strouhal_number(Re) for Re in Re_range]

ax1.semilogx(Re_range, St_range, 'k-', linewidth=3, label='Empirical Correlation')

# Mark test cases
for case in test_cases:
    ax1.plot(case['Re'], case['St'], 'o', color=case['color'],
            markersize=15, markeredgewidth=2, markeredgecolor='black',
            label=f"{case['regime']} (Re={case['Re']:.1e}, St={case['St']:.3f})")

# Add regime shading
ax1.axvspan(50, 300, alpha=0.15, color='blue')
ax1.axvspan(300, 2e5, alpha=0.15, color='green')
ax1.axvspan(2e5, 3.5e6, alpha=0.15, color='orange')
ax1.axvspan(3.5e6, 1e7, alpha=0.15, color='red')

# Regime labels
ax1.text(150, 0.32, 'Laminar', fontsize=11, ha='center', weight='bold')
ax1.text(5e3, 0.32, 'Subcritical', fontsize=11, ha='center', weight='bold')
ax1.text(7e5, 0.32, 'Critical/\nSupercritical', fontsize=10, ha='center', weight='bold')
ax1.text(5e6, 0.32, 'Transcritical', fontsize=11, ha='center', weight='bold')

ax1.set_xlabel('Reynolds Number (Re)', fontsize=13, weight='bold')
ax1.set_ylabel('Strouhal Number (St)', fontsize=13, weight='bold')
ax1.set_title('Vortex Shedding Frequency vs Reynolds Number', fontsize=14, weight='bold')
ax1.grid(True, alpha=0.3, which='both')
ax1.legend(fontsize=11, loc='lower right')
ax1.set_xlim([50, 1e7])
ax1.set_ylim([0, 0.35])

# ============================================================
# Plot 2: Turbulent Eddy Viscosity Coefficient
# ============================================================
ax2 = fig.add_subplot(gs[1, 0])

Re_range_fine = np.logspace(2, 7, 200)
C_eddy_range = [compute_eddy_viscosity_coeff(Re) for Re in Re_range_fine]

ax2.semilogx(Re_range_fine, C_eddy_range, 'purple', linewidth=3)

# Mark test cases
for case in test_cases:
    C_eddy = compute_eddy_viscosity_coeff(case['Re'])
    ax2.plot(case['Re'], C_eddy, 'o', color=case['color'],
            markersize=12, markeredgewidth=2, markeredgecolor='black')
    ax2.annotate(f"{case['regime']}\nC_eddy={C_eddy:.3f}",
                xy=(case['Re'], C_eddy),
                xytext=(20, 10), textcoords='offset points',
                fontsize=9, ha='left',
                bbox=dict(boxstyle='round', facecolor=case['color'], alpha=0.3),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))

# Threshold line
ax2.axvline(1000, color='red', linestyle='--', linewidth=2, alpha=0.7,
           label='Turbulence Activation (Re=1000)')

ax2.set_xlabel('Reynolds Number (Re)', fontsize=12, weight='bold')
ax2.set_ylabel('Eddy Viscosity Coefficient $C_{eddy}$', fontsize=12, weight='bold')
ax2.set_title('Turbulent Eddy Viscosity vs Re', fontsize=13, weight='bold')
ax2.grid(True, alpha=0.3, which='both')
ax2.legend(fontsize=10)
ax2.set_xlim([50, 1e7])
ax2.set_ylim([0, 0.12])

# ============================================================
# Plot 3: Shedding Period Comparison
# ============================================================
ax3 = fig.add_subplot(gs[1, 1])

Re_test = [100, 500, 1000, 5000, 1e4, 1e5, 1e6, 6.7e6]
St_test = [get_strouhal_number(Re) for Re in Re_test]
T_shed = [D / (St * U_inf) for St in St_test]
f_shed = [1.0 / T for T in T_shed]

colors = ['blue' if Re <= 1000 else 'red' for Re in Re_test]

ax3.bar(range(len(Re_test)), f_shed, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)

# Add values on bars
for i, (Re, f) in enumerate(zip(Re_test, f_shed)):
    ax3.text(i, f + 0.005, f'{f:.3f} Hz', ha='center', fontsize=9, weight='bold')

ax3.set_xticks(range(len(Re_test)))
ax3.set_xticklabels([f'{Re:.0e}' for Re in Re_test], rotation=45, ha='right')
ax3.set_xlabel('Reynolds Number (Re)', fontsize=12, weight='bold')
ax3.set_ylabel('Shedding Frequency [Hz]', fontsize=12, weight='bold')
ax3.set_title('Vortex Shedding Frequency Across Regimes', fontsize=13, weight='bold')
ax3.grid(True, alpha=0.3, axis='y')
ax3.set_ylim([0, max(f_shed) * 1.2])

# Add legend
blue_patch = Rectangle((0, 0), 1, 1, fc='blue', alpha=0.7, edgecolor='black')
red_patch = Rectangle((0, 0), 1, 1, fc='red', alpha=0.7, edgecolor='black')
ax3.legend([blue_patch, red_patch], ['Laminar (Re<1000)', 'Turbulent (Re>1000)'],
          fontsize=10, loc='upper left')

# ============================================================
# Plot 4: Feature Activation Summary
# ============================================================
ax4 = fig.add_subplot(gs[2, :])

features = ['Eddy\nViscosity', 'Stochastic\nShedding', 'Core\nSaturation',
           'Re-Dependent\nStrouhal']
feature_labels = ['Eddy Viscosity', 'Stochastic Shedding', 'Core Saturation', 'Re-Dependent St']

# Status for each test case
status_laminar = [0, 0, 0, 1]  # Only dynamic St is active
status_turbulent = [1, 1, 1, 1]  # All features active

x = np.arange(len(features))
width = 0.35

bars1 = ax4.bar(x - width/2, status_laminar, width, label='Laminar (Re=100)',
               color='blue', alpha=0.7, edgecolor='black', linewidth=1.5)
bars2 = ax4.bar(x + width/2, status_turbulent, width, label='Transcritical (Re=6.7M)',
               color='red', alpha=0.7, edgecolor='black', linewidth=1.5)

# Add status labels
for i, (v1, v2) in enumerate(zip(status_laminar, status_turbulent)):
    if v1 == 1:
        ax4.text(i - width/2, v1 + 0.05, 'ON', ha='center', fontsize=10, weight='bold')
    else:
        ax4.text(i - width/2, 0.5, 'OFF', ha='center', fontsize=10, weight='bold', color='gray')

    if v2 == 1:
        ax4.text(i + width/2, v2 + 0.05, 'ON', ha='center', fontsize=10, weight='bold')

ax4.set_ylabel('Feature Status', fontsize=12, weight='bold')
ax4.set_xlabel('Turbulence Feature', fontsize=12, weight='bold')
ax4.set_title('Turbulence Feature Activation: Laminar vs Transcritical Regime',
             fontsize=13, weight='bold')
ax4.set_xticks(x)
ax4.set_xticklabels(features, fontsize=11)
ax4.set_yticks([0, 1])
ax4.set_yticklabels(['OFF', 'ON'], fontsize=11, weight='bold')
ax4.legend(fontsize=11, loc='upper right')
ax4.grid(True, alpha=0.3, axis='y')
ax4.set_ylim([0, 1.3])

# Add text box with key differences
textstr = '\n'.join([
    'Key Differences:',
    '─' * 50,
    f'• Laminar (Re=100):',
    f'  - St = 0.167  →  f = 0.051 Hz  →  T = 19.5 s',
    f'  - Molecular viscosity only (ν = {U_inf*D/100:.2e} m²/s)',
    f'  - Deterministic shedding',
    f'  - Coherent vortices persist 20-30D downstream',
    '',
    f'• Transcritical (Re=6.7M):',
    f'  - St = 0.300  →  f = 0.092 Hz  →  T = 10.9 s',
    f'  - Eddy viscosity >> molecular (ν_t ≈ 0.1·U·D)',
    f'  - Stochastic shedding (±10-20% fluctuations)',
    f'  - Wake spreads 5-10× faster, vortex breakdown',
])

props = dict(boxstyle='round', facecolor='wheat', alpha=0.8, edgecolor='black', linewidth=2)
fig.text(0.02, 0.02, textstr, fontsize=10, verticalalignment='bottom',
        bbox=props, family='monospace')

plt.suptitle('Laminar vs Turbulent Flow: Comparison of Wake Dynamics',
            fontsize=16, weight='bold', y=0.995)

plt.tight_layout(rect=[0, 0.15, 1, 0.99])
plt.savefig('comparison_laminar_vs_turbulent.png', dpi=300, bbox_inches='tight')
print("\n[OK] Saved comparison_laminar_vs_turbulent.png")
print("\nComparison plot shows:")
print("  1. Strouhal number variation across Re regimes")
print("  2. Eddy viscosity coefficient scaling")
print("  3. Shedding frequency at test Re values")
print("  4. Feature activation status (laminar vs turbulent)")
print("\nKey findings:")
print(f"  • Laminar (Re=100):   f_shed = 0.051 Hz  (T = 19.5 s)")
print(f"  • Turbulent (Re=6.7M): f_shed = 0.092 Hz  (T = 10.9 s)")
print(f"  • Frequency increased by 80% due to higher St at transcritical Re")
