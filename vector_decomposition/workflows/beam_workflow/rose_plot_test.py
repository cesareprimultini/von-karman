import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from windrose import WindroseAxes

# =================CONFIGURATION=================
# Replace with your actual CSV file path
CSV_FILE_PATH = r'vector_decomposition\workflows\beam_workflow\P2_data.csv'

# Replace with the exact column names in your CSV for speed and direction
# Direction should be in degrees (0-360), Speed in m/s
SPEED_COL = 'uc'
DIR_COL = 'current_dir'

# Define the specific speed bins for the legend coloring.
# Adjust these values to match the ranges seen in the reference image if needed.
SPEED_BINS = [0, 0.2, 0.4, 0.6, 0.8, 1.0]

# Choose a colormap that matches the style (Spectral_r goes from blue to red)
COLOR_MAP = cm.Spectral_r
# ===============================================


def plot_cumulative_rose(csv_path, speed_col, dir_col, bins, cmap):
    try:
        # 1. Load Data
        df = pd.read_csv(csv_path)

        # Basic validation
        if speed_col not in df.columns or dir_col not in df.columns:
            print(f"Error: Columns '{speed_col}' and/or '{dir_col}' not found in CSV.")
            return

        # Extract data and drop any NaN values just in case
        data = df[[speed_col, dir_col]].dropna()
        speed = data[speed_col].values
        direction = data[dir_col].values

        if len(speed) == 0:
             print("Error: No valid data found after dropping NaNs.")
             return

        # 2. Setup Figure and Windrose Axes
        fig = plt.figure(figsize=(10, 10))
        rect = [0.1, 0.1, 0.8, 0.8]
        # We use WindroseAxes specifically for this type of polar plot
        ax = WindroseAxes(fig, rect)
        fig.add_axes(ax)

        # 3. Create the stacked bar plot (Rose Plot)
        # normed=True ensures the radial axis is in Percentage (%)
        # opening=0.9 determines the width of the sectors (1.0 is touching, <1 leaves gaps)
        # edgecolor='white' adds the thin lines between stacked segments
        ax.bar(direction, speed,
               normed=True,
               opening=0.9,
               bins=bins,
               nsector=12,
               cmap=cmap,
               edgecolor='white',
               linewidth=0.5)

        # 4. Customize Axes and Grid Style
        # Setting the radial grid lines (percentages).
        # Find a suitable max percentage based on data, or hardcode (e.g., up to 30%)
        # Here we set ticks every 5% up to 30%
        yticks = np.arange(5, 51, 5)
        ax.set_yticks(yticks)
        ax.set_yticklabels([f'{i}%' for i in yticks])

        # Set 12 angular grid lines at 30° intervals
        angles = np.arange(0, 360, 30)
        labels = ['N', '30', '60', 'E', '120', '150', 'S', '210', '240', 'W', '300', '330']
        ax.set_thetagrids(angles, labels=labels)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)

        # Style grid lines
        ax.grid(True, linestyle='--', alpha=0.7, color='grey')

        # Style labels (N, E, S, W and percentages)
        ax.tick_params(axis='x', pad=15, labelsize=12,labelcolor='black', colors='black') # Direction labels
        ax.tick_params(axis='y', labelsize=10, labelcolor='grey', colors='grey') # Percentage labels


        # 5. Add Legend
        # We place it outside the plot area to the bottom right
        legend = ax.set_legend(title="$U_c$ Speed [m/s]",
                               loc='lower right',
                               bbox_to_anchor=(1.2, 0.05),
                               fancybox=True, shadow=False)

        # Adjust legend text size
        plt.setp(legend.get_texts(), fontsize=10)
        plt.setp(legend.get_title(), fontsize=11, fontweight='bold')

        # Add a title
        plt.title("Cumulative Current Rose Plot", y=1.08, fontsize=14, fontweight='bold')

        # 6. Save output
        output_filename = "current_rose_plot_cumulative.png"
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Plot generated successfully and saved to: {output_filename}")
        # plt.show() # Uncomment if you want to display it immediately

    except FileNotFoundError:
        print(f"Error: File not found at {csv_path}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Run the function
if __name__ == "__main__":
    # Prerequisite check: Ensure windrose library is installed.
    # If missing, run: pip install windrose
    try:
        import windrose
    except ImportError:
        print("Error: This script requires the 'windrose' library.")
        print("Please install it using: pip install windrose")
        exit()

    plot_cumulative_rose(CSV_FILE_PATH, SPEED_COL, DIR_COL, SPEED_BINS, COLOR_MAP)