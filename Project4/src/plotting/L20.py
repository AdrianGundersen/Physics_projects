import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
data = np.loadtxt(ROOT / "data/outputs/L20_random.txt")  # add delimiter=',' if CSV
x = np.arange(len(data))

plt.plot(x, data); plt.xlabel("Index"); plt.ylabel("Value"); plt.title("L20 random"); plt.grid(True)
out = ROOT / "data/plots/L20_random.pdf"
out.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out, format="pdf", bbox_inches="tight"); print(f"Saved: {out}")
plt.show()