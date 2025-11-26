import numpy as np

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

            parts = [p for p in line.split(",") if p != ""]
            if len(parts) != 2:
                raise ValueError(f"Line does not have two comma-separated values: {line}")
            re = float(parts[0])
            im = float(parts[1])
            psi2 = (re*re + im*im)
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

    M = prob_fields[0].shape 
    h = 1.0/(M[0] - 1)
    prob_fields = np.array(prob_fields) * h * h  # normalize

    return prob_fields