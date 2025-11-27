import numpy as np
import matplotlib.pyplot as plt
import os
from io_python import read_prob_file

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

output_dir = "output/figures"
os.makedirs(output_dir, exist_ok=True)


def plot_screen_distribution(prob_fields, dt, T, x_screen=0.8, L=1.0, slits="single"):
    """
    Plots p(y | x=x_screen, t=T).

    """
    fileneme = f"detection_probability_{slits}_slit.pdf"

    t_index = int(T / dt)
    field = prob_fields[t_index]      
    M = field.shape[0]
    dx = L / (M - 1)

    # index to find x_screen
    j_screen = int(round(x_screen / dx))

    # probability along the screen
    line_prob = field[j_screen, :]     # |psi(x_screen, y_i, T)|^2

    # p(y | x_screen, T)
    p_y_given_x = line_prob / (np.sum(line_prob))
    print("Sum p(y|x) =", np.sum(p_y_given_x))

    y = (np.arange(M)) * dx


    plt.figure()
    plt.plot(y, p_y_given_x)
    plt.xlabel(r"$y$")
    plt.ylabel(r"$p(y \mid x=0.8;\, t=0.002)$")
    plt.tight_layout()
    out_path = os.path.join(output_dir, fileneme)
    plt.savefig(out_path)
    plt.close()
    print("Saved figure:", out_path)


if __name__ == "__main__":
    filename = "output/wavefunction_p8_200.txt"  # juster sti
    prob_fields = read_prob_file(filename)

    slits = "200grid"

    dt = 2.5e-5
    T = 0.002
    L = 1.0
    x_screen = 0.8

    plot_screen_distribution(prob_fields, dt, T, x_screen=x_screen, L=L, slits=slits)
