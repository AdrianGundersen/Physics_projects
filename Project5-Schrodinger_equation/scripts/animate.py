import numpy as np
import matplotlib
matplotlib.use("Agg") 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from io_python import read_prob_file

# plotting style
plt.rcParams.update({
    "text.usetex": False, 
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "cm",        
    'font.size': 16,
    'figure.figsize': (6, 4),
    'axes.titlesize': 20,
    'axes.labelsize': 20,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'lines.linewidth': 2.0,
    'legend.fontsize': 18,
    'figure.dpi': 300,
})

def plot_timestep(prob_fields, t_index=0):
    field = prob_fields[t_index]
    # field = np.sqrt(field) # to visualize better
    plt.figure()
    im = plt.imshow(field, origin="lower")
    plt.colorbar(im, label=r"$|\psi|^2$")
    plt.title(f"Probability density, timestep {t_index}")
    plt.xlabel("j")
    plt.ylabel("i")
    plt.tight_layout()
    plt.show()


def animate_prob(prob_fields, dt=1.0, frame_stride=1):
    """
    Animate the probability density fields.

    prob_fields : list of 2D arrays, in time order.
    dt          : physical time step between stored frames.
    frame_stride: use every N-th frame (1 = all frames).
    """
    fontsize = 12

    # Time indices used in the animation
    frame_indices = list(range(0, len(prob_fields), frame_stride))

    # Global colour normalization across all timesteps
    vmax = max(np.max(field) for field in prob_fields)
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=vmax)

    # Create figure/axes
    fig, ax = plt.subplots()

    # Initial frame
    first_field = prob_fields[frame_indices[0]]
    img = ax.imshow(first_field, origin="lower", cmap="viridis", norm=norm)

    # Axis labels
    ax.set_xlabel(r"$y$", fontsize=fontsize)
    ax.set_ylabel(r"$x$", fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)

    # Colourbar
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label(r"$|v_{{ij}}^{{n}}|^2$", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    fig.tight_layout()

    # Time text (in axes coordinates)
    time_txt = ax.text(
        0.95,
        0.95,
        f"t = {0.0:.3e}",
        color="white",
        ha="right",
        va="top",
        fontsize=fontsize,
        transform=ax.transAxes,
    )

    def update(frame_idx):
        field = prob_fields[frame_idx]
        img.set_data(field)

        t = frame_idx * dt
        time_txt.set_text(f"t = {t:.3e}")

        return img, time_txt

    anim = FuncAnimation(
        fig,
        update,
        frames=frame_indices,
        interval=50,   # ms between frames
        repeat=False,
        blit=False,
    )

    return anim


if __name__ == "__main__":
    filename_suffix = "wavefunction" # problem
    filename_interfix = "" # or double, multi etc.
    filename = f"output/{filename_suffix}_{filename_interfix}.txt"  # adjust path if needed
    prob_fields = read_prob_file(filename)

    # Example: static plot of a single timestep
    # plot_timestep(prob_fields, t_index=40)

    # Animate the full time evolution
    # Set dt to your actual simulation time step if you want a physical time axis
    anim = animate_prob(prob_fields, dt=2.5e-5, frame_stride=1)

    # To save instead of (or in addition to) showing:
    # anim.save("probability_animation.mp4", writer="ffmpeg", bitrate=10000, fps=15)
    anim.save(f"probability_animation_{filename_interfix}.gif", writer="pillow", fps=15)
    plt.show()

