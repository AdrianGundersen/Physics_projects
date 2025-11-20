# Project 4

## Description
Simulation of the 2D Ising model using the Metropolis algorithm.

## Physical Model
We consider a 2D square $L \times L$ lattice of spins $s_i = \pm 1$ with periodic boundary conditions and without an external magnetic field. The total energy of the system is given by:
$$
E = -J \sum_{\langle i,j \rangle}^N s_i s_j, \qquad J>0,\; N = L^2
$$
where the sum runs over all nearest-neighbor pairs \(\langle i,j \rangle\). The magnetization is defined as:
$$
M = \sum_{i=1}^N s_i
$$

We also have the specific heat capacity $c_V$ and magnetic susceptibility $\chi_s$ defined as:

$$
c_V = \frac{1}{N k_B T^2} \left( \langle E^2 \rangle - \langle E \rangle^2 \right)
$$
$$
\chi_s = \frac{1}{N k_B T} \left( \langle M^2 \rangle - \langle |M| \rangle^2 \right)
$$
We will however work in units where \(J = 1\) and \(k_B = 1\), so these expressions simplify to and only look at per spin quantitites.
We will use this to study phase transitions in the system as we vary the temperature \(T\). 

## Project Structure
```text
Project4/
├── Makefile                 # Build targets for simulators/plotting helpers
├── README.md                # Project overview and quickstart
├── PROJECT_PLAN.md          # Internal scope/timeline notes
├── configs/                 # JSON configs for reproducible runs
│   ├── 2x2_test.json        # Analytical validation
│   ├── L20.json             # Burn-in diagnostics
│   ├── multiple_walkers.json # Scaling test (cores x walkers)
│   └── t_sweep.json         # Temperature sweep definition
├── data/                    # Generated artifacts (git-ignored)
│   ├── figures/             # Saved plots
│   └── outputs/             # Time series, QR lattices, JSON aggregates
├── include/                 # C++ headers
│   ├── omp_rng.hpp
│   └── ising/
│       ├── lattice.hpp
│       ├── metropolis.hpp
│       ├── mcmc_run.hpp
│       ├── model.hpp
│       ├── observables.hpp
│       └── io/
│           ├── json_util.hpp
│           └── write_files.hpp
├── src/                     # C++ implementations
│   ├── omp_rng.cpp
│   └── ising/
│       ├── lattice.cpp
│       ├── metropolis.cpp
│       ├── mcmc_run.cpp
│       ├── observables.cpp
│       └── io/
│           └── write_files.cpp
├── apps/                    # Entry-point executables
│   ├── Ln.cpp               # Production simulator (single T or sweeps)
│   └── runtime.cpp          # Timing harness for scaling studies
├── plotting/                # Python analysis + visualization
│   ├── 2x2.py               # Numerical vs analytical L=2 comparison
│   ├── L20.py               # Burn-in inspection for L=20 data
│   ├── exp_vals.py          # CLI helper for averaged observables/errors
│   ├── plot_cv_chi.py       # Cv/chi curves and Tc extrapolation
│   ├── plot_speedup.py      # Parallel speedup plot from measurements
│   └── qr_code.py           # Convert saved lattices to QR-style images
├── bin/                     # Compiled executables (generated)
│   ├── Ln
│   └── runtime
└── test/
    └── JSON_example/        # Reference output / fixtures
```

## Plot Config overview

All other parameters like `measure_sweeps: 1` and `seed: 67` are fixed for all simulations. `"burn_in_sweeps": 10000` unless said otherwise

| Plot/Script & App                  | Config file(s)             | Lattice size(s) | Temperature setup                              | Walkers x threads                 | Monte Carlo cycles (`total_sweeps`) | Notes / TODO |
|-----------------------------------|----------------------------|-----------------|-----------------------------------------------|-----------------------------------|-------------------------------------|--------------|
| `plotting/L20.py` + `Ln.cpp`      | `configs/L20.json`         | 20              | Single-T runs at 1.0 and 2.4                   | 1 x 1                             | 1e6                  | Burn-in traces + histograms |
| `plotting/2x2.py` + `Ln.cpp`      | `configs/2x2_test.json`    | 2               | Both `use_Trange=true` and `false`, \(T\in[1,4]\) | 1 x 1                             | 1e6 and 1e7                         | Compare numerics vs analytical solution |
| `plotting/plot_cv_chi.py` + `Ln.cpp` | `configs/t_sweep.json`   | 5-120           | `use_Trange=true`, see `data/tsweep_finale.json` | 12 x 12                           | see `data/tsweep_finale.json`                | Extract Cv/chi peaks + Tc estimate |
| `plotting/plot_speedup.py` + `apps/runtime.cpp` | `configs/multiple_walkers.json` | 20 | Single T                                      | 48 walkers, threads \([1,16])     | 1e4, no burn-in.                                 | Update runtimes/core counts manually |
| `plotting/qr_code.py` + `Ln.cpp`  | `configs/multiple_walkers.json` w/ `qr=true` | 40 | Single T = 2.377                               | 12 x 12                           | 1e6               | Produces QR-style lattice snapshot |



## Installation and Dependencies


- **C++20 compiler** (e.g. g++ 11 or higher)
- **Armadillo** 
- **OpenMP** 
- **Python 3.10+**
- **nhlomann/json** (header-only C++ JSON library) - included in    `include/` folder
  - `numpy`
  - `matplotlib`


### Ubuntu/Debian
```bash
sudo apt update
sudo apt install g++ make libarmadillo-dev python3 python3-pip
pip3 install numpy matplotlib
```

Other systems not tested, but should work with similar steps.


## Usage & Compilation
To see how make targets are defined, check the `Makefile` or run:
```bash
make help
```

### Using JSON config files
The executables in the `apps/` folder read parameters from JSON config files located in the `configs/` folder. You can create your own config files by copying and modifying the existing templates. The argument `use_Trange` tells the program to run for a single temeprature of not. 

**Example**
How we did the first search for peaks at $L=20$.
```json
{
    "model": // Model parameters
    {
        "J": 1.0, // Interaction strength
        "double_count" : false, // Whether to double count energy
        "spin_config": "random" // Initial spin configuration: "all_up", "all_down", "random"
    },
    "lattice": // Lattice parameters
    {
        "L" : 20 // Lattice size (LxL)
    },
    "simulation": // Simulation parameters (Metropolis)
    {
        "seed": 67, // Mother seed for random number generator
        "temperature": 2.0, // Temperature [J/kB], using this temperature if "use_Trange": false
        "use_Trange": true, // Bool too run over multiple temeratures (makes temperature parameter obsolete)
        "Trange": {
            "Tmin": 2.1,    // Starting temerature
            "Tmax": 2.4,    // End temeprature
            "Tsteps": 11    // Number of temerature steps from Tmin to Tmax
        },
        "total_steps": "N", // Total Monte Carlo steps ("N" means number af spins but can also be any integer)
        "burn_in_sweeps": 10000, // Burn-in-cycles before starting mearuring
        "measure_sweeps": 1, // Measure observables every n sweeps
        "total_sweeps": 1e6, // Total number of sweeps after burn-in
        "cores": 12, // Number of CPU cores to use
        "walkers": 12, // Number of independent walkers (should be = n * cores for best performance)
        "qr": false // Bool to write lattice spins to file, only works if "use_Trange": false.
    },
    "write_to_file": 
    {
        "enabled": true, // whether to write results to file
        "observables": ["Cv", "chi"], // what observables to write [e, m] for txt [cv, chi] for json
        "type": "json", // file type: "txt", "json" (only txt & json supported currently)
        "delimiter": ",", // delimiter for txt files
        "precision": 6, // number of decimal places
        "output_filename": "default", // output filename (if "default", uses autogenerated name)
        "header": false, // only for txt
        "average_or_concatenate": "concatenate" // average takes average over all wakers, concatenate merges all data from walkers into one file (average NEVER recommended as it loses info, but was used for testing purposes)
    }
}
```

**Caveats**
- For writing to json:
  - Only T_range mode is supported. 
  - Will not care about "observables" field. Always writes (Cv, chi, avg_eps, avg_mabs, avg_eps2, avg_mabs2, avg_eps3, avg_mabs3, avg_eps4, avg_mabs4).

- For writing to txt:
    - T_range mode not supported.
    - Only eps and m supported. Gives energy and magnetization per spin.
    - Recommended to use only one walker for ordered data.

### Compile
To compile the code to `bin/`, use:
```bash
make all
```

or to compile only specific:
```bash
make <target>
```

### Run
To run the code locally, use:
```bash 
make run
```
and the cinfig file `multiple_walkers.json` we be automatically used.

**Raspberry Pi**
We have also used a raspberry pi 5 for extra computational power. To run raspberry pi, use:
```bash
make PI_HOST=<your_pi_host>
```

or set the `PI_HOST` variable in the `Makefile` and run:
```bash
make run-pi
```
Depending on your version and OS, you might need to change `PI_ARCH` variable in the `Makefile` to match your architecture.


### Plotting

### What different scripts do

**Simulation apps**
- `apps/Ln.cpp`: Main entry point for Ising simulations; reads a JSON config, runs either a single-temperature Metropolis sweep or a full temperature range, writes observables, and optionally stores a QR-style lattice snapshot.
- `apps/runtime.cpp`: Minimal harness that runs one simulation from a JSON config and prints elapsed time together with the requested OpenMP core count—used for scaling studies and sanity checks.

**Plotting and analysis utilities (Python)**
- `plotting/L20.py`: Inspects burn-in for \(L=20\) runs by plotting raw energies, running means, and post burn-in histograms for different initial spin configurations.
- `plotting/plot_cv_chi.py`: Aggregates JSON sweep outputs for many lattice sizes, plots \(C_V(T)\) and \(\chi(T)\), and extracts peak locations/uncertainties used to estimate the critical temperature.
- `plotting/2x2.py`: Compares \(L=2\) simulations against the analytical solution; tracks running relative errors, observable-vs-temperature curves, and saves all diagnostic figures.
- `plotting/exp_vals.py`: CLI helper (`python exp_vals.py <file> <L> <T>`) that reads a two-column data file, computes averaged observables, and reports relative errors vs the exact \(L=2\) solution.
- `plotting/plot_speedup.py`: Uses hard-coded timing measurements to plot measured vs ideal parallel speedup and prints the corresponding time-saving factor.
- `plotting/qr_code.py`: Turns a saved lattice snapshot (written when `qr=true`) into a black/white QR-style image for visualizing spin domains.

## Important parameters for reproducability


## Future work
- Generalize Metropolis for any model.
- Make a better write_to_file script that can have more options and raises errors when unsupported options are chosen. Also structure it into multiple functions for better readability.

## Declaration of Use of Generative AI

In this scientific work, generative artificial intelligence (AI) has been used. All data and personal information have been processed in accordance with the University of Oslo’s regulations, and we, as the authors of the document, take full responsibility for its content, claims, and references. An overview of the use of generative AI is provided below.

**Summary**

- **Tool(s) used:** [OpenAI ChatGPT (GPT-5)](https://chatgpt.com/)
- **Use:**
  - Generating boilerplate code (e.g., plotting with `matplotlib`).
  - Generating Makefiles.
  - Formatting `README.md`.
  - Checking language for clarity and grammar; general proofreading.
  - Generating brief code documentation, comments, and variable names to improve readability and adhere to conventions.
  - Creating tables in proper TeX format.
  - Brainstorming ideas for code optimization and efficiency.
