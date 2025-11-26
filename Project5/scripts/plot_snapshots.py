import numpy as np
import matplotlib.pyplot as plt

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
    field = np.abs(psi) ** 2 / psi.size

    plt.figure()
    im = plt.imshow(field, origin="lower", extent=[0, field.shape[1]*dx, 0, field.shape[0]*dx])
    plt.colorbar(im, label=r"$|\psi_{ij}|^2$")
    plt.xlabel("y")
    plt.ylabel("x")
    plt.tight_layout()
    plt.show()


def plot_re_im_timestep(psi_fields, t_index=0, dt=1.0, dx=1.0):
    """
    Plot colourmaps of Re(psi_ij) and Im(psi_ij) for a given timestep index.
    """
    psi = psi_fields[t_index]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    im_re = axes[0].imshow(psi.real, origin="lower", extent=[0, psi.shape[1]*dx, 0, psi.shape[0]*dx])
    plt.colorbar(im_re, ax=axes[0], label=r"$\Re(u_{ij})$")
    axes[0].set_xlabel("y")
    axes[0].set_ylabel("x")

    im_im = axes[1].imshow(psi.imag, origin="lower", extent=[0, psi.shape[1]*dx, 0, psi.shape[0]*dx])
    plt.colorbar(im_im, ax=axes[1], label=r"$\Im(u_{ij})$")
    axes[1].set_xlabel("y")
    axes[1].set_ylabel("x")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    filename = "output/wavefunction.txt"  # adjust path if needed
    psi_fields = read_wavefile(filename)
    
    dt = 2.5e-5  # time step size in seconds
    T = 0
    t_index = int(T / dt)

    L = 1.0 # box size
    M = psi_fields[0].shape[0]
    dx = L / M

    plot_prob_timestep(psi_fields, t_index=t_index, dt=dt, dx=dx)
    plot_re_im_timestep(psi_fields, t_index=t_index, dt=dt, dx=dx)