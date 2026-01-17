"""Visualization functions for von Kármán vortex simulation results."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from von_karman_simulator import compute_velocity_field, apply_transforms


def plot_velocity_field(results_df, timestep_index, cylinders, flow_params, plot_config):
    """
    Reconstruct and plot total velocity field (freestream + vortices).

    Parameters
    ----------
    results_df : pd.DataFrame
        DataFrame from VortexAmp.run()
    timestep_index : int
        Row index to plot
    cylinders : list[dict]
        Cylinder geometry: [{'x': x, 'y': y, 'D': D}, ...]
    flow_params : dict
        {'rotation_angle': float, 'flow_angle': float}
    plot_config : dict
        {'x_range': (min, max), 'y_range': (min, max), 'grid_size': int,
         'arrow_skip': int, 'dpi': int, 'filename': str}
    """
    row = results_df.iloc[timestep_index]
    vortex_field = row['vortex_field']

    if vortex_field is None:
        print(f"No vortex field saved at timestep {timestep_index}. Choose a timestep where save_interval triggered.")
        return

    time = row['time']
    U_inf = row['U_inf']
    Re = row['Re']

    print(f"Computing velocity field for plot at t={time:.1f}s...")

    x_range = plot_config.get('x_range', (-80, 150))
    y_range = plot_config.get('y_range', (-50, 50))
    grid_size = plot_config.get('grid_size', 600)
    arrow_skip = plot_config.get('arrow_skip', 4)
    dpi = plot_config.get('dpi', 600)
    filename = plot_config.get('filename', 'velocity_field.png')

    x = np.linspace(x_range[0], x_range[1], grid_size)
    y = np.linspace(y_range[0], y_range[1], grid_size)
    X, Y = np.meshgrid(x, y)

    U, V = compute_velocity_field(X, Y, vortex_field, U_inf, flow_params['flow_angle'])
    vel_mag = np.sqrt(U**2 + V**2)

    fig, ax = plt.subplots(figsize=(16, 8), dpi=300)

    contour = ax.contourf(X, Y, vel_mag, levels=100, cmap='viridis')
    plt.colorbar(contour, ax=ax, label='Velocity Magnitude [m/s]')

    skip = arrow_skip
    scale_reference = np.percentile(vel_mag, 95)
    scale_dynamic = scale_reference * 80
    ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
              U[::skip, ::skip], V[::skip, ::skip],
              color='white', alpha=0.7, scale=scale_dynamic, width=0.0005)

    for cyl in cylinders:
        x_cyl, y_cyl = apply_transforms(cyl['x'], cyl['y'],
                                        flow_params['rotation_angle'],
                                        flow_params['flow_angle'])
        circle = Circle((x_cyl, y_cyl), cyl['D']/2,
                        color='gray', fill=True, zorder=10,
                        edgecolor='black', linewidth=2)
        ax.add_patch(circle)

    for v in vortex_field:
        color = 'ro' if v['gamma'] > 0 else 'bo'
        ax.plot(v['x'], v['y'], color, markersize=2, alpha=0.4)

    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    D_ref = max(cyl['D'] for cyl in cylinders)
    ax.set_title(f'Vortex Wake (D={D_ref}m, Re={Re:.2e}, t={time:.1f}s)')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    print(f"Saved {filename}")
    plt.close()


def plot_vortex_perturbations(results_df, timestep_index, cylinders, flow_params, plot_config):
    """
    Reconstruct and plot vortex-induced perturbations ONLY (no freestream).

    This computes velocity field with U_inf=0, showing pure vortex effects.
    Useful for isolating vortex contributions before adding cylinder potential
    flow or other amplification factors.

    Parameters
    ----------
    results_df : pd.DataFrame
        DataFrame from VortexAmp.run()
    timestep_index : int
        Row index to plot
    cylinders : list[dict]
        Cylinder geometry: [{'x': x, 'y': y, 'D': D}, ...]
    flow_params : dict
        {'rotation_angle': float, 'flow_angle': float}
    plot_config : dict
        {'x_range': (min, max), 'y_range': (min, max), 'grid_size': int,
         'arrow_skip': int, 'dpi': int, 'filename': str}
    """
    row = results_df.iloc[timestep_index]
    vortex_field = row['vortex_field']

    if vortex_field is None:
        print(f"No vortex field saved at timestep {timestep_index}. Choose a timestep where save_interval triggered.")
        return

    time = row['time']
    Re = row['Re']

    print(f"Computing vortex perturbations for plot at t={time:.1f}s...")

    x_range = plot_config.get('x_range', (-80, 150))
    y_range = plot_config.get('y_range', (-50, 50))
    grid_size = plot_config.get('grid_size', 600)
    arrow_skip = plot_config.get('arrow_skip', 4)
    dpi = plot_config.get('dpi', 600)
    filename = plot_config.get('filename', 'vortex_perturbations.png')

    x = np.linspace(x_range[0], x_range[1], grid_size)
    y = np.linspace(y_range[0], y_range[1], grid_size)
    X, Y = np.meshgrid(x, y)

    U, V = compute_velocity_field(X, Y, vortex_field, U_inf=0.0, flow_angle=0.0)
    vel_mag = np.sqrt(U**2 + V**2)

    fig, ax = plt.subplots(figsize=(16, 8), dpi=300)

    contour = ax.contourf(X, Y, vel_mag, levels=100, cmap='viridis')
    plt.colorbar(contour, ax=ax, label='Perturbation Velocity Magnitude [m/s]')

    skip = arrow_skip
    scale_reference = np.percentile(vel_mag, 95)
    scale_dynamic = scale_reference * 80
    ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
              U[::skip, ::skip], V[::skip, ::skip],
              color='white', alpha=0.7, scale=scale_dynamic, width=0.0005)

    for cyl in cylinders:
        x_cyl, y_cyl = apply_transforms(cyl['x'], cyl['y'],
                                        flow_params['rotation_angle'],
                                        flow_params['flow_angle'])
        circle = Circle((x_cyl, y_cyl), cyl['D']/2,
                        color='gray', fill=True, zorder=10,
                        edgecolor='black', linewidth=2)
        ax.add_patch(circle)

    for v in vortex_field:
        color = 'ro' if v['gamma'] > 0 else 'bo'
        ax.plot(v['x'], v['y'], color, markersize=2, alpha=0.4)

    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    D_ref = max(cyl['D'] for cyl in cylinders)
    ax.set_title(f'Vortex Perturbations (D={D_ref}m, Re={Re:.2e}, t={time:.1f}s)')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    print(f"Saved {filename}")
    plt.close()


def plot_velocity_history(results_df, probe_indices=None, mark_shedding=False, t_start=None, t_end=None):
    """
    Plot velocity time history for specified probes.

    Parameters
    ----------
    results_df : pd.DataFrame
        DataFrame from simulation
    probe_indices : list[int]
        List of probe indices to plot (None = all)
    mark_shedding : bool
        Whether to mark shedding events (vertical lines when vortices shed)
    t_start : float, optional
        Start time for plot range (None = use data minimum)
    t_end : float, optional
        End time for plot range (None = use data maximum)
    """
    if probe_indices is None:
        probe_cols = [col for col in results_df.columns if col.startswith('probe_') and col.endswith('_ux')]
        probe_indices = [int(col.split('_')[1]) for col in probe_cols]

    if len(probe_indices) == 0:
        print("No measurement points found in results.")
        return

    # Filter dataframe by time range
    time_mask = np.ones(len(results_df), dtype=bool)
    if t_start is not None:
        time_mask &= (results_df['time'] >= t_start)
    if t_end is not None:
        time_mask &= (results_df['time'] <= t_end)

    results_filtered = results_df[time_mask]

    if len(results_filtered) == 0:
        print(f"No data found in time range [{t_start}, {t_end}].")
        return

    print("Generating velocity measurement plot...")
    fig, ax = plt.subplots(figsize=(10, 6))

    time = results_filtered['time'].values

    for idx in probe_indices:
        ux = results_filtered[f'probe_{idx}_ux'].values
        uy = results_filtered[f'probe_{idx}_uy'].values
        vmag = results_filtered[f'probe_{idx}_vmag'].values

        ax.plot(time, ux, label=f'Probe {idx}: $u_x$', linewidth=1.0, alpha=0.7)
        ax.plot(time, uy, label=f'Probe {idx}: $u_y$', linewidth=1.0, alpha=0.7)
        ax.plot(time, vmag, label=f'Probe {idx}: $|v|$', linewidth=1.5)

    if mark_shedding:
        n_vortices = results_filtered['n_vortices'].values
        shed_indices = np.where(np.diff(n_vortices) > 0)[0]
        shed_times = time[shed_indices]
        for t_shed in shed_times:
            ax.axvline(t_shed, color='red', alpha=0.15, linewidth=0.5, linestyle='-')

    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Velocity [m/s]', fontsize=12)
    ax.set_title('Velocity Time History at Measurement Points', fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, ncol=2)
    ax.set_xlim([time[0], time[-1]])

    plt.tight_layout()
    plt.savefig('velocity_history.png', dpi=300, bbox_inches='tight')
    print("Saved velocity_history.png")
    plt.close()


def plot_freestream_history(results_df):
    """
    Plot U_inf and Re vs time (for time-varying cases).

    Parameters
    ----------
    results_df : pd.DataFrame
        DataFrame from simulation
    """
    time = results_df['time'].values
    U_inf = results_df['U_inf'].values
    Re = results_df['Re'].values

    if np.allclose(U_inf, U_inf[0]):
        print("Freestream velocity is constant. Skipping freestream history plot.")
        return

    print("Generating freestream velocity plot...")
    fig, ax = plt.subplots(figsize=(10, 4))

    ax.plot(time, U_inf, 'b-', linewidth=2, label='$U_\\infty(t)$')
    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Freestream Velocity [m/s]', fontsize=12, color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.set_title('Time-Varying Freestream Velocity', fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([time[0], time[-1]])

    ax2 = ax.twinx()
    ax2.plot(time, Re, 'r--', linewidth=2, label='Re(t)')
    ax2.set_ylabel('Reynolds Number', fontsize=12, color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    ax2.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=11)

    plt.tight_layout()
    plt.savefig('freestream_history.png', dpi=300, bbox_inches='tight')
    print("Saved freestream_history.png")
    plt.close()


def create_animation(results_df, cylinders, flow_params, plot_config,
                     output_file='wake_animation.mp4', frame_skip=1,
                     sample_frames=30, cleanup_frames=True, fps=30):
    """
    Create animated visualization using optimized frame pre-rendering.

    This function samples frames to determine consistent colorbar limits, pre-renders
    frames to PNG files, then compiles with ffmpeg for fast video generation.

    Parameters
    ----------
    results_df : pd.DataFrame
        DataFrame from simulation
    cylinders : list[dict]
        Cylinder geometry: [{'x': x, 'y': y, 'D': D}, ...]
    flow_params : dict
        {'rotation_angle': float, 'flow_angle': float}
    plot_config : dict
        {'x_range': tuple, 'y_range': tuple, 'grid_size': int, 'arrow_skip': int}
    output_file : str
        Output video filename (use .mp4, .avi, or .gif extension)
    frame_skip : int
        Use every Nth snapshot (1 = use all, 2 = use every other, etc.)
    sample_frames : int
        Number of frames to sample for determining colorbar limits (default: 30)
    cleanup_frames : bool
        Remove temporary PNG frame files after creating video (default: True)
    fps : int
        Frames per second for output video (default: 30)

    Returns
    -------
    None

    Notes
    -----
    Requires ffmpeg to be installed and available in system PATH.
    For GIF output, the pillow writer will be used instead of ffmpeg.
    """
    import os
    import shutil

    snapshot_indices = results_df[results_df['vortex_field'].notna()].index.tolist()

    if len(snapshot_indices) == 0:
        print("No vortex field snapshots found. Cannot create animation.")
        return

    snapshot_indices = snapshot_indices[::frame_skip]
    total_frames = len(snapshot_indices)

    print(f"Creating animation with {total_frames} frames (skip={frame_skip})...")

    x_range = plot_config.get('x_range', (-80, 150))
    y_range = plot_config.get('y_range', (-50, 50))
    grid_size = plot_config.get('grid_size', 400)
    arrow_skip = plot_config.get('arrow_skip', 6)

    x = np.linspace(x_range[0], x_range[1], grid_size)
    y = np.linspace(y_range[0], y_range[1], grid_size)
    X, Y = np.meshgrid(x, y)

    sample_step = max(1, total_frames // sample_frames)
    sample_indices = snapshot_indices[::sample_step]

    print(f"Sampling {len(sample_indices)} frames to determine colorbar limits...")
    vel_mag_min = np.inf
    vel_mag_max = -np.inf

    from tqdm import tqdm

    for idx in tqdm(sample_indices, desc="Sampling frames", unit="frame"):
        row = results_df.iloc[idx]
        vortex_field = row['vortex_field']
        U_inf = row['U_inf']
        U, V = compute_velocity_field(X, Y, vortex_field, U_inf, flow_params['flow_angle'])
        vel_mag = np.sqrt(U**2 + V**2)
        vel_mag_min = min(vel_mag_min, vel_mag.min())
        vel_mag_max = max(vel_mag_max, vel_mag.max())

    print(f"Sampled velocity range: [{vel_mag_min:.2f}, {vel_mag_max:.2f}] m/s")

    frames_dir = 'temp_animation_frames'
    os.makedirs(frames_dir, exist_ok=True)

    for frame_idx, idx in enumerate(tqdm(snapshot_indices, desc="Rendering frames", unit="frame")):
        row = results_df.iloc[idx]
        vortex_field = row['vortex_field']
        time = row['time']
        U_inf = row['U_inf']
        Re = row['Re']

        fig, ax = plt.subplots(figsize=(16, 8), dpi=100)

        U, V = compute_velocity_field(X, Y, vortex_field, U_inf, flow_params['flow_angle'])
        vel_mag = np.sqrt(U**2 + V**2)

        contour = ax.contourf(X, Y, vel_mag, levels=50, cmap='viridis',
                             vmin=vel_mag_min, vmax=vel_mag_max)

        cbar = plt.colorbar(contour, ax=ax, label='Velocity Magnitude [m/s]')

        skip = arrow_skip
        scale_reference = np.percentile(vel_mag, 95)
        scale_dynamic = scale_reference * 80
        ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
                  U[::skip, ::skip], V[::skip, ::skip],
                  color='white', alpha=0.7, scale=scale_dynamic, width=0.0005)

        for cyl in cylinders:
            x_cyl, y_cyl = apply_transforms(cyl['x'], cyl['y'],
                                            flow_params['rotation_angle'],
                                            flow_params['flow_angle'])
            circle = Circle((x_cyl, y_cyl), cyl['D']/2,
                            color='gray', fill=True, zorder=10,
                            edgecolor='black', linewidth=2)
            ax.add_patch(circle)

        for v in vortex_field:
            color = 'ro' if v['gamma'] > 0 else 'bo'
            ax.plot(v['x'], v['y'], color, markersize=2, alpha=0.4)

        ax.set_xlim(x_range)
        ax.set_ylim(y_range)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        D_ref = max(cyl['D'] for cyl in cylinders)
        ax.set_title(f'Vortex Wake (D={D_ref}m, Re={Re:.2e}, t={time:.1f}s)')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)

        frame_file = os.path.join(frames_dir, f'frame_{frame_idx:06d}.png')
        plt.savefig(frame_file, dpi=100, bbox_inches='tight')
        plt.close(fig)

    print(f"\nAll frames rendered. Compiling video with ffmpeg...")

    file_ext = os.path.splitext(output_file)[1].lower()

    if file_ext == '.gif':
        print("Creating GIF with pillow (this may take a moment)...")
        import matplotlib.animation as animation

        fig = plt.figure(figsize=(16, 8))
        frames = []
        for frame_idx in range(total_frames):
            frame_file = os.path.join(frames_dir, f'frame_{frame_idx:06d}.png')
            img = plt.imread(frame_file)
            im = plt.imshow(img, animated=True)
            plt.axis('off')
            frames.append([im])

        ani = animation.ArtistAnimation(fig, frames, interval=1000/fps, blit=True, repeat_delay=1000)
        ani.save(output_file, writer='pillow')
        plt.close()
    else:
        try:
            import subprocess
            frame_pattern = os.path.join(frames_dir, 'frame_%06d.png')

            cmd = [
                'ffmpeg',
                '-y',
                '-framerate', str(fps),
                '-i', frame_pattern,
                '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2',
                '-c:v', 'libx264',
                '-pix_fmt', 'yuv420p',
                '-crf', '18',
                output_file
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode != 0:
                print(f"ffmpeg error: {result.stderr}")
                print("Falling back to pillow writer...")
                import matplotlib.animation as animation
                fallback_file = os.path.splitext(output_file)[0] + '.gif'
                print(f"Saving as {fallback_file} instead...")
                fig = plt.figure(figsize=(16, 8))
                frames = []
                for frame_idx in range(total_frames):
                    frame_file = os.path.join(frames_dir, f'frame_{frame_idx:06d}.png')
                    img = plt.imread(frame_file)
                    im = plt.imshow(img, animated=True)
                    plt.axis('off')
                    frames.append([im])
                ani = animation.ArtistAnimation(fig, frames, interval=1000/fps, blit=True)
                ani.save(fallback_file, writer='pillow')
                plt.close()
                output_file = fallback_file
        except FileNotFoundError:
            print("ffmpeg not found. Falling back to pillow writer...")
            import matplotlib.animation as animation
            fallback_file = os.path.splitext(output_file)[0] + '.gif'
            print(f"Saving as {fallback_file} instead...")
            fig = plt.figure(figsize=(16, 8))
            frames = []
            for frame_idx in range(total_frames):
                frame_file = os.path.join(frames_dir, f'frame_{frame_idx:06d}.png')
                img = plt.imread(frame_file)
                im = plt.imshow(img, animated=True)
                plt.axis('off')
                frames.append([im])
            ani = animation.ArtistAnimation(fig, frames, interval=1000/fps, blit=True)
            ani.save(fallback_file, writer='pillow')
            plt.close()
            output_file = fallback_file

    if cleanup_frames:
        print(f"Cleaning up temporary frames...")
        shutil.rmtree(frames_dir)
    else:
        print(f"Temporary frames saved in {frames_dir}/")

    print(f"Animation saved to {output_file}")
    print(f"  Total frames: {total_frames}")
    print(f"  Frame rate: {fps} fps")
    print(f"  Duration: {total_frames/fps:.1f} seconds")
