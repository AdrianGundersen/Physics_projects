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


def plot(prob_fields):
    """
    Plot total prob density against time steps
    """
    t_steps = len(prob_fields)
    total_probs = []
    for t in range(t_steps):
        total_prob = np.sum(prob_fields[t])
        total_probs.append(total_prob)
    plt.figure()
    plt.plot(range(t_steps), total_probs, marker='o')
    plt.xlabel("Timestep")
    plt.ylabel("Total Probability")
    plt.title("Total Probability vs Timestep")
    plt.grid()
    plt.tight_layout()
    plt.show()
    return total_probs

if __name__ == "__main__":
    filename = "output/wavefunction.txt" 
    prob_fields = read_prob_file(filename)
    total_probs = plot(prob_fields)
    print("Total probabilities at each timestep:", total_probs)