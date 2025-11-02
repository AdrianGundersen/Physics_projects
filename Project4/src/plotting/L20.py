import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Finds project root
ROOT = Path(__file__).resolve().parents[2]

# Spin configs and temps
spin_configs = ["all_up", "all_down", "random"]
temperatures = [1.0, 2.40]  # flere kan legges til


fig_dir = ROOT / "data/figures/L20"
fig_dir.mkdir(parents=True, exist_ok=True)

for spin_config in spin_configs:
    for T in temperatures:
        # data/outputs/L20_T=1.000000_spin=all_up.txt
        data_path = ROOT / f"data/outputs/L20_T={T:0.6f}_spin={spin_config}.txt"

        if not data_path.exists():
            print(f"Warning: {data_path} does not exist.")
            continue

        try:
            data = np.loadtxt(data_path)  # legg til delimiter=',' hvis du senere skriver CSV
        except Exception as e:
            print(f"Error: Could not read {data_path}: {e}")
            continue

        x = np.arange(len(data))

        # Title and labels
        title = f"L20 {spin_config} T={T:0.6f}"
        xlabel = "Sweep index"
        ylabel = rf"$\epsilon$"

        # --- Plotting ---
        plt.figure()
        plt.plot(x, data)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(True)
        out_pdf = fig_dir / f"L20_{spin_config}_T{T:0.6f}.pdf"
        plt.savefig(out_pdf, format="pdf", bbox_inches="tight")
        plt.close()
        print(f"Saved: {out_pdf}")