# Project 3

## Description
Simulations for a Penning trap. The code models motion of charged particles in an electromagnetic field. There are three main simulations:

- **Single-particle dynamics** — Simulates motion with **Forward Euler** and **RK4** integrators, compared against analytical solutions.
- **Few-particle interactions** — Models Coulomb interaction between two charged particles in the trap.
- **Many-particle parameter sweeps** — Investigates trapping stability by sweeping over driving frequency and amplitude.
- **Floquet analysis** — Determines regions of stability based on the Mathieu equation formulation. (incomplete)
- **Parallelization (OpenMP)** — Multi-core acceleration for systems with many particles.



## Project structure
```text
Project3/
├─ include/
│  ├─ constants.hpp
│  ├─ integrator.hpp
│  ├─ parameters.hpp
│  ├─ Particle.hpp
│  └─ PenningTrap.hpp
├─ scripts/
│  ├─ animate.py
│  ├─ floquet_table.py
│  ├─ plot_pos_vel.py
│  ├─ plot_single.py
│  └─ plot_trapped_vs_omega.py
├─ src/
│  ├─ integrator.cpp
│  ├─ Particle.cpp
│  ├─ PenningTrap.cpp
│  └─ few_particles.cpp
├─ N_particles.cpp
├─ single_particle.cpp
├─ Makefile
└─ README.md
```

## Installation and Dependencies


- **C++20 compiler** (e.g. g++ 11 or higher)
- **Armadillo** (version 12 or higher)
- **OpenMP** (can be disabled)
- **Python 3.10+**
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

Change parameters in 
```bash
include/parameters.hpp 
```
as needed. The parameters are well documented in the file.




### Compile

Build everything
```bash
make
```

Build in debug mode
```bash
make debug
```

Build with sanitizers (useful for debugging memory issues)
```bash
make sanitize
```

Run individual simulations
```bash
make run-s
make run-few
make run-n
```

Clean objects and executables:
```bash
make clean
```
### Run
```bash
From Project root:
./build/bin/single_particle
./build/bin/few_particles
./build/bin/N_particles
```

Or build, then run:
```bash
make run
```

Display help:
```bash
make help
```

Toggle OpenMP (off)
```bash
make OMP=0
```


### What different scripts do:
- **Example** -

### Plotting
Have to run the simulations first to generate data files. Then run the python scripts in the scripts/ folder to generate plots and animations. Example:
```bash
python3 scripts/plot_single.py
```
Requirements to replicate plots in the report:
- `plot_single.py` run `make run-s` for step sizes $N=4000, 8000, 16000, 32000$, with total time $50\mu s$.
- `plot_pos_vel.py` run `make run-few` for step size $N=40{,}000$, total time $50\mu s$ with and without coulumb interactions. 
- `plot_trapped_vs_omega.py` first run `make run-few` for step size $N=40{,}000$, total time $500\mu s$ without coulumb interactions Change the window for $\omega_V$ to `w_min = 0.2, w_min = 0.2, w_max = 2.5, w_step = 0.001`. Then in `plot_trapped_vs_omega.py` set enter the variables for the window $\omega_V$. To run narrow windows eneter window parameters end set `w_step = 0.0005` and run `make run-few` with and without coulumb forces. 



Remember also to change filenames and change the parameters in the scripts to match those in parameters.hpp.


| Script | Description | Input | Output |
|---------|--------------|--------|---------|
| `plot_single.py` | Compares analytical vs numerical trajectories for one particle | Data files from `single_particle.exe` | Time evolution & relative error plots |
| `plot_pos_vel.py` | Plots positions and velocities of multiple particles | Few/N-particle data files | Phase space and trajectory plots |
| `plot_trapped_vs_omega.py` | Plots trapped fraction vs driving frequency | Parameter sweep results | Stability diagrams |
| `floquet_table.py` | Generates a Floquet/Mathieu stability table | Parameters inside | LaTeX-formatted stability tables |
| `animate.py` | Creates animations of particle motion | Positional data files | `.gif` or `.mp4` animations |


## Future work
- Add unit tests
- Add Verlet/Leapfrog integrator
- Further evolve the flouquet_table.py script to generate a proper stability diagram
- Make a parameters.json file to read parameters from instead of hardcoding them in parameters.hpp. This would also make it easier to keep parameters in sync between the C++ code and the python scripts.
