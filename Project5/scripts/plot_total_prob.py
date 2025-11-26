import numpy as np
import matplotlib.pyplot as plt
import os

# plotting style
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 12,
    'figure.figsize': (6, 4),
    'axes.titlesize': 14,
    'axes.labelsize': 20,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'lines.linewidth': 2.0,
    'legend.fontsize': 14,
    'figure.dpi': 300,
})


# make output/figures directory
output_dir = "output/figures"
os.makedirs(output_dir, exist_ok=True)






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
    L = 1.0
    h = L/(grid_size - 1)
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
            psi2 = (re*re + im*im) * (h * h)  # normalize here
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


def plot(prob_fields, output_dir_prefix=None):
    """
    Plot total prob density against time steps
    """
    t_steps = len(prob_fields)
    total_probs = []
    for t in range(t_steps):
        total_prob = np.sum(prob_fields[t])
        total_probs.append(total_prob)

    p0 = total_probs[0]
    dp =  np.array(total_probs) -p0
    plt.figure()
    plt.plot(range(t_steps), dp, marker='o')
    plt.xlabel("Timestep")
    plt.ylabel("Total Probability")
    plt.title("Total Probability vs Timestep")
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"{output_dir_prefix}_total_prob.png")
    plt.show()

    # histogram
    plt.figure()
    plt.hist(dp, bins=40)
    plt.xlabel(r"$P(t)-P(0)$")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(f"{output_dir_prefix}_hist.png")
    plt.show()
    return total_probs

if __name__ == "__main__":
    filename_with_slit = "output/wavefunction_slit_p7.txt"
    filename_no_slit = "output/wavefunction_no_slit.txt"

    output_prefix_with_slit = output_dir + "/total_prob_with_slit"
    output_prefix_no_slit = output_dir + "/total_prob_no_slit"

    prob_fields_with_slit = read_prob_file(filename_with_slit)
    total_probs_with_slit = plot(prob_fields_with_slit, output_prefix_with_slit)
    prob_fields_no_slit = read_prob_file(filename_no_slit)
    total_probs_no_slit = plot(prob_fields_no_slit, output_prefix_no_slit)