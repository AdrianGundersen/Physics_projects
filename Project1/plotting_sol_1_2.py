import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('diff_eq_sol.txt')

x = data[:, 0]
u = data[:, 1]

plt.plot(x, u)
plt.show()