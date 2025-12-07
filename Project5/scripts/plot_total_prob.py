import numpy as np
import matplotlib.pyplot as plt
import os
from io_python import read_prob_file

# plotting style
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 14,
    'figure.figsize': (6, 4),
    'axes.titlesize': 16,
    'axes.labelsize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'lines.linewidth': 2.0,
    'legend.fontsize': 18,
    'figure.dpi': 300,
})

# make output/figures directory
output_dir = "output/figures"
os.makedirs(output_dir, exist_ok=True)

def plot(prob_fields, output_dir_prefix=None):
    """
    Plot total prob density against time steps
    """
    t_steps = len(prob_fields)
    total_probs = []
    for t in range(t_steps):
        total_prob = np.sum(prob_fields[t])
        total_probs.append(total_prob)

    p0_exp = 1.0
    p0 = total_probs[0]
    dp =  np.array(total_probs) - p0_exp
    plt.figure()
    plt.plot(range(t_steps), dp, marker='o')
    plt.xlabel("Timestep")
    plt.ylabel("Total Probability")
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"{output_dir_prefix}_total_prob.pdf")
    #plt.show()

    # histogram
    mean = np.mean(dp)
    std = np.std(dp)
    print(f"For file {output_dir_prefix}: P0 = {p0}, mean = {mean}, std = {std}")

    plt.figure()
    plt.hist(dp, bins=40)
    plt.axvline(mean, color='r', linestyle='--')
    plt.axvline(mean + std, color='g', linestyle=':')
    plt.axvline(mean - std, color='g', linestyle=':')
    plt.xlabel(r"$P(t)-P(0)$")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(f"{output_dir_prefix}_hist.pdf")
    #plt.show()
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