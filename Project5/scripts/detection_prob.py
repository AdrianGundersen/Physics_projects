import numpy as np
import matplotlib.pyplot as plt
import os
from io_python import read_prob_file

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
    dx = L / M

    # index to find x_screen
    j_screen = int(x_screen / dx)

    # probability along the screen
    line_prob = field[j_screen, :]     # |psi(x_screen, y_i, T)|^2

    # p(y | x_screen, T)
    p_y_given_x = line_prob / np.sum(line_prob)

    y = (np.arange(M) + 0.5) * dx


    plt.figure()
    plt.plot(y, p_y_given_x)
    plt.xlabel(r"$y$")
    plt.ylabel(r"$p(y \mid x=0.8;\, t=0.002)$")
    plt.title("Deteksjonssannsynlighet langs skjermen ved $x=0.8$")
    plt.tight_layout()
    out_path = os.path.join(output_dir, fileneme)
    plt.savefig(out_path)
    plt.close()
    print("Saved figure:", out_path)


if __name__ == "__main__":
    filename = "output/wavefunction_p8.txt"  # juster sti
    prob_fields = read_prob_file(filename)

    slits = "double"

    dt = 2.5e-5
    T = 0.002
    L = 1.0
    x_screen = 0.8

    plot_screen_distribution(prob_fields, dt, T, x_screen=x_screen, L=L, slits=slits)
