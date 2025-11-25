import numpy as np
import matplotlib.pyplot as plt

def read_prob_file(filename):
    """
    Read file with blocks of the form

    Timestep 0:
    Re0, Im0,
    Re_1, Im_1,
    ...

    separated by blank lines.
    Returns: list of 2D numpy arrays [prob_t0, prob_t1, ...]
    """
    blocks = []
    current_vals = []
    grid_size = 200  # adjust if needed

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

            #read values
            parts = line.split(",")
            if len(parts) != 2:
                raise ValueError(f"Line does not have two comma-separated values: {line}")
            re = float(parts[0])
            im = float(parts[1])
            psi2 = (re*re + im*im) / (grid_size*grid_size)  # normalize here
            current_vals.append(psi2)

    if current_vals:
        blocks.append(current_vals)

    prob_fields = []
    for vals in blocks:
        n = len(vals)
        M = int(np.sqrt(n))
        if M * M != n:
            raise ValueError(f"Block has {n} values, which is not a perfect square.")
        arr = np.array(vals).reshape(M, M)
        prob_fields.append(arr)

    return prob_fields

def plot_timestep(prob_fields, t_index=0):
    """
    Plot probability density for a given timestep index.
    """
    field = prob_fields[t_index]

    plt.figure()
    im = plt.imshow(field, origin="lower")
    plt.colorbar(im, label=r"$|\psi|^2$")
    plt.title(f"Probability density, timestep {t_index}")
    plt.xlabel("j")
    plt.ylabel("i")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    filename = "output/wavefunction.txt"  # adjust path if needed
    prob_fields = read_prob_file(filename)
    # Example: plot timestep 0
    plot_timestep(prob_fields, t_index=40)
    # Example: plot the last timestep
    # plot_timestep(prob_fields, t_index=len(prob_fields) - 1)