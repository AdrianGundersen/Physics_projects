import matplotlib.pyplot as plt
import numpy as np
import os

# Path to this script's directory
script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir, "diff_eq_sol.txt")

data = np.loadtxt(file_path)

x = data[:, 0]
u = data[:, 1]

plt.plot(x, u)
plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Solution to the Differential Equation")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("Plot2")
plt.show()
plt.close()