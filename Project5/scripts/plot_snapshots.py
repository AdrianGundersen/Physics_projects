import numpy as np
import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt
import os

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

# make output/figures directory
output_dir = "output/figures"
os.makedirs(output_dir, exist_ok=True)


def read_wavefile(filename):
    """
    Read file with blocks of the form

    Timestep 0:
    Re0, Im0,
    Re1, Im1,
    ...

    separated by blank lines.

    Returns: list of 2D numpy arrays of complex numbers [psi_t0, psi_t1, ...]
    """
    blocks = []
    current_vals = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                if current_vals:
                    blocks.append(current_vals)
                    current_vals = []
                continue

            if line.startswith("Timestep"):
                continue

            parts = [p for p in line.split(",") if p != ""]
            if len(parts) != 2:
                raise ValueError(f"Line does not have two comma-separated values: {line}")

            re = float(parts[0])
            im = float(parts[1])
            current_vals.append(re + 1j * im)

    if current_vals:
        blocks.append(current_vals)

    psi_fields = []
    for vals in blocks:
        n = len(vals)
        M = int(np.sqrt(n))
        if M * M != n:
            raise ValueError(f"Block has {n} values, which is not a perfect square.")
        arr = np.array(vals, dtype=np.complex128).reshape(M, M)
        psi_fields.append(arr)

    return psi_fields


def plot_prob_timestep(psi_fields, t_index=0, dt=1.0, dx=1.0):
    """
    Plot probability density |psi|^2 for a given timestep index.
    Normalised by the total number of grid points.
    """
    psi = psi_fields[t_index]
    field = np.abs(psi) ** 2
    plt.figure()
    im = plt.imshow(field, origin="lower", extent=[0, field.shape[1]*dx, 0, field.shape[0]*dx])
    plt.colorbar(im, label=rf"$|v_{{ij}}^{{{t_index}}}|^2$")
    plt.xlabel(r"$y$")
    plt.ylabel(r"$x$")
    plt.tight_layout()
    plt.savefig(output_dir + f"/probability_density_timestep={t_index}.pdf")
    print("Saved figure:", output_dir + f"/probability_density_timestep={t_index}.pdf")
    # plt.show()

def plot_prob_timestep_zoom_square(psi_fields, t_index=0, dt=1.0, dx=1.0, zoom=False, preffix=""):
    """
    Plot probability density |psi|^2 for a given timestep index.
    Normalised by the total number of grid points.
    """
    psi_fields = np.sqrt(psi_fields)
    psi = psi_fields[t_index]
    field = np.abs(psi) ** 2
    plt.figure()
    im = plt.imshow(field, origin="lower", extent=[0, field.shape[1]*dx, 0, field.shape[0]*dx])
    if zoom==True:
        plt.xlim(0.3, 0.7)
        plt.ylim(0.3, 0.7) 
    plt.colorbar(im, label=rf"$|v_{{ij}}^{{{t_index}}}|$")
    plt.xlabel(r"$y$")
    plt.ylabel(r"$x$")
    plt.tight_layout()
    plt.savefig(output_dir + f"/probability_density_timestep_{preffix}={t_index}.pdf")
    print("Saved figure:", output_dir + f"/probability_density_timestep_{preffix}={t_index}.pdf")
    # plt.show()


def plot_re_im_timestep(psi_fields, t_index=0, dt=1.0, dx=1.0):
    """
    Plot colourmaps of Re(psi_ij) and Im(psi_ij) for a given timestep index.
    """
    psi_fields = np.sqrt(psi_fields)  # rescale for better visibility
    psi = psi_fields[t_index]

    # 2 rows, 1 column
    fig, axes = plt.subplots(2, 1, figsize = (6,8), sharex=True, sharey=True)

    im_re = axes[0].imshow(
        psi.real,
        origin="lower",
        extent=[0, psi.shape[1] * dx, 0, psi.shape[0] * dx]
    )
    plt.colorbar(im_re, ax=axes[0], label=rf"$\sqrt{{\Re(v_{{ij}}^{{{t_index}}})}}$")
    axes[0].set_ylabel(r"$x$")
    axes[0].set_title(rf"Real part $\,\sqrt{{\Re(v_{{ij}}^{{{t_index}}})}}$")

    im_im = axes[1].imshow(
        psi.imag,
        origin="lower",
        extent=[0, psi.shape[1] * dx, 0, psi.shape[0] * dx]
    )
    plt.colorbar(im_im, ax=axes[1], label=rf"$\sqrt{{\Im(v_{{ij}}^{{{t_index}}})}}$")
    axes[1].set_xlabel(r"$y$")
    axes[1].set_ylabel(r"$x$")
    axes[1].set_title(rf"Imaginary part $\,\sqrt{{\Im(v_{{ij}}^{{{t_index}}})}}$")

    plt.tight_layout()
    plt.savefig(output_dir + f"/re_im_timestep={t_index}.pdf")
    print("Saved figure:", output_dir + f"/re_im_timestep={t_index}.pdf")
    # plt.show()

def plot_re_timestep(psi_fields, t_index=0, dt=1.0, dx=1.0):
    """
    Plot colourmap of Re(psi_ij) zoomed in around slit for a given timestep index.
    """
    psi_fields = np.sqrt(psi_fields)  # rescale for better visibility
    psi = psi_fields[t_index]

    fig, ax = plt.subplots(figsize=(6, 4))
    im_re = ax.imshow(
        psi.real,
        origin="lower",
        extent=[0, psi.shape[1] * dx, 0, psi.shape[0] * dx],
        interpolation="bilinear",   
        aspect="equal",
    )

    ax.set_ylim(0.2, 0.8)
    ax.set_xlim(0.2, 0.8)

    fig.colorbar(im_re, ax=ax, label=rf"$\sqrt{{\Re(v_{{ij}}^{{{t_index}}})}}$")
    ax.set_xlabel(r"$y$")
    ax.set_ylabel(r"$x$")
    ax.set_title(rf"Real part $\,\sqrt{{\Re(v_{{ij}}^{{{t_index}}})}}$")

    fig.tight_layout()
    fig.savefig(output_dir + f"/re_zoom_timestep={t_index}.pdf")
    print("Saved figure:", output_dir + f"/re_zoom_timestep={t_index}.pdf")
    plt.close(fig)




if __name__ == "__main__":
    filename = "output/wavefunction_2slit.txt"  # adjust path if needed  
    filename_single = "output/wavefunction_1slit.txt"  # adjust path if needed  
    psi_fields = read_wavefile(filename)
    psi_fields_single = read_wavefile(filename_single)
    
    dt = 2.5e-5  # time step size in seconds
    T = np.array([0.0, 0.001, 0.002])

    t_index_list = (T / dt).astype(int)

    L = 1.0 # box size
    M = psi_fields[0].shape[0]
    dx = L / M

    zoom = False
    for t_index in t_index_list:
        if t_index == t_index_list[1]:
            plot_re_timestep(psi_fields_single, t_index=t_index, dt=dt, dx=dx)
            zoom = True
        else:
            zoom = False
        plot_prob_timestep_zoom_square(psi_fields, t_index=t_index, dt=dt, dx=dx, preffix="2slit")
        plot_re_im_timestep(psi_fields, t_index=t_index, dt=dt, dx=dx)
        plot_prob_timestep_zoom_square(psi_fields_single, t_index=t_index, dt=dt, dx=dx, zoom=zoom, preffix="1slit")