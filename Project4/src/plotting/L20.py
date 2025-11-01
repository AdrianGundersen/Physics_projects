import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Load data (1D list of values)
data = np.loadtxt('/home/adriang/FYS3150/Project4/data/outputs/L20_random.txt')

# x = index (0..N-1), y = data values
x = np.arange(len(data), dtype=int)
y = data

# Make plot
plt.figure()
plt.plot(x, y)
plt.xlabel("Index")
plt.ylabel("Value")
plt.title("L20 random")
plt.grid(True)

# Save as PDF (ensure directory exists)
out_dir = Path('/home/adriang/FYS3150/Project4/data/figures')
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / 'L20_random.pdf'
plt.savefig(out_path, format='pdf', bbox_inches='tight')

plt.show()
print(f"Saved plot to: {out_path}")