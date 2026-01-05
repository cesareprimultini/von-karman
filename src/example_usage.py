"""Example usage of the modular von Kármán vortex simulation."""

from von_karman_simulator import VortexAmp
import plotting


# EXAMPLE: Basic constant velocity simulation

print("=" * 60)
print("EXAMPLE 1: Constant velocity simulation")
print("=" * 60)

# Define geometry
cylinders = [
    {'x': -20, 'y': +20, 'D': 4.88},
    {'x': +20, 'y': +20, 'D': 4.88},
    {'x': -20, 'y': -20, 'D': 4.88},
    {'x': +20, 'y': -20, 'D': 4.88},
]

# Create simulator with 6 measurement points
sim = VortexAmp(
    cylinders=cylinders,
    nu=4.88e-6,
    flow_angle_metocean=270.0,
    rotation_angle=30.0,
    dt=0.01,
    measurement_points=[(32.0, 5.0)],
    save_interval=5.0  # Save vortex field every 5 seconds (can be None for no saving)
)

# Run with constant velocity
results = sim.run(
    velocity_mode='constant',
    U_inf=1.5,
    total_time=200.0  # Shorter for example
)

# Save results
sim.save_results('test_01.pkl')
print("\n")


# Load and verify results

print("=" * 60)
print("EXAMPLE: Load and verify results")
print("=" * 60)

# Load results
results = VortexAmp.load_results('test_01.pkl')

print(f"Results DataFrame shape: {results.shape}")
print(f"Columns: {list(results.columns)}")
print(f"Time range: {results['time'].min():.1f}s - {results['time'].max():.1f}s")
print(f"Number of snapshots with vortex field: {results['vortex_field'].notna().sum()}")
print("\n")


# EXAMPLE 3: Plot velocity field at specific time

print("=" * 60)
print("EXAMPLE 3: Plot velocity field")
print("=" * 60)

# Find a timestep where vortex field was saved
snapshot_indices = results[results['vortex_field'].notna()].index.tolist()
if len(snapshot_indices) > 0:
    # Plot at middle of simulation
    plot_idx = snapshot_indices[len(snapshot_indices) // 2]

    # Plot total velocity field (freestream + vortices)
    plotting.plot_velocity_field(
        results,
        timestep_index=plot_idx,
        cylinders=cylinders,
        flow_params={'rotation_angle': 30.0, 'flow_angle': 0.0},
        plot_config={
            'x_range': (-80, 150),
            'y_range': (-50, 50),
            'grid_size': 600,
            'arrow_skip': 4,
            'dpi': 600,
            'filename': 'example_velocity_field.png'
        }
    )

    # Plot vortex perturbations only (no freestream)
    plotting.plot_vortex_perturbations(
        results,
        timestep_index=plot_idx,
        cylinders=cylinders,
        flow_params={'rotation_angle': 30.0, 'flow_angle': 0.0},
        plot_config={
            'x_range': (-80, 150),
            'y_range': (-50, 50),
            'grid_size': 600,
            'arrow_skip': 4,
            'dpi': 600,
            'filename': 'example_vortex_perturbations.png'
        }
    )
print("\n")


# EXAMPLE 4: Plot velocity time histories

print("=" * 60)
print("EXAMPLE 4: Plot velocity time histories")
print("=" * 60)

# Plot all probes
plotting.plot_velocity_history(results, probe_indices=None, mark_shedding=True)
print("\n")


# EXAMPLE 5: Time-varying velocity and direction

print("=" * 60)
print("EXAMPLE 5: Time-varying velocity and direction")
print("=" * 60)

def ramp_velocity(t):
    """Linear velocity ramp"""
    return 1.0 + 0.01 * t

def oscillating_direction(t):
    """Oscillating flow direction (metocean convention)"""
    import numpy as np
    return 270.0 + 30.0 * np.sin(2 * np.pi * t / 25.0)

sim2 = VortexAmp(
    cylinders=cylinders,
    nu=4.88e-6,
    flow_angle_metocean=270.0,
    rotation_angle=30.0,
    dt=0.01,
    measurement_points=[(32.0, 5.0)],
    save_interval=10.0
)

results2 = sim2.run(
    velocity_mode='function',
    velocity_function=ramp_velocity,
    flow_angle_mode='function',
    flow_angle_function=oscillating_direction,
    total_time=50.0
)

plotting.plot_freestream_history(results2)
print("\n")


# EXAMPLE 6: Custom post-processing

print("=" * 60)
print("EXAMPLE 6: Custom post-processing")
print("=" * 60)

# Extract data from specific time
t_target = 25.0
row = results[results['time'].round(2) == t_target].iloc[0]

print(f"At t={t_target}s:")
print(f"  U_inf = {row['U_inf']:.3f} m/s")
print(f"  Re = {row['Re']:.2e}")
print(f"  St = {row['St']:.3f}")
print(f"  Number of vortices = {row['n_vortices']}")
print(f"  Probe 0 velocity = ({row['probe_0_ux']:.3f}, {row['probe_0_uy']:.3f}) m/s")

# Analyze vortex shedding frequency
import numpy as np
from scipy import signal

# Use probe 0 time history
time = results['time'].values
vmag = results['probe_0_vmag'].values

# Skip initial transient
t_start = 10.0
idx_start = np.argmin(np.abs(time - t_start))

# Compute FFT
dt = time[1] - time[0]
frequencies, psd = signal.welch(vmag[idx_start:], fs=1/dt, nperseg=1024)

# Find dominant frequency
idx_peak = np.argmax(psd)
f_dominant = frequencies[idx_peak]
St_measured = f_dominant * cylinders[0]['D'] / row['U_inf']

print(f"\nShedding frequency analysis:")
print(f"  Dominant frequency = {f_dominant:.3f} Hz")
print(f"  Measured Strouhal = {St_measured:.3f}")
print(f"  Theoretical Strouhal = {row['St']:.3f}")
print("\n")


print("=" * 60)
print("All examples complete!")
print("=" * 60)
