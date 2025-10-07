# Project 3

## Description
Simulations for a Penning trap. The code models motion of charged particles in an electromagnetic field. There are three main simulations:

- **Single particle dynamics** with both Forward Euler and RK4 integrators and comparison to analytical result.
- **Few-particle interactions** with simulation of two particles in the trap.
- **Parameter sweeps (many particles)** with driving frequency and amplitude added.



### Project structure
```text
Project3/
├─ include/
│  ├─ constants.hpp
│  ├─ integrator.hpp
│  ├─ parameters.hpp
│  ├─ Particle.hpp
│  └─ PenningTrap.hpp
├─ script/
│  ├─ animate.py
│  ├─ plot_pos_vel.py
│  └─plot_single.py
├─ src/
│  ├─ integrator.cpp
│  ├─ Particle.cpp
│  └─ PenningTrap.cpp
├─ few_particles.cpp
├─ N_particles.cpp
├─ single_particle.cpp
├─ build/
├─ Makefile
└─ README.md
```

## Dependencies

- **C++20 compiler**
- **Armadillo**   
- **OpenMP**
- **Python 3.10+**
  - `numpy`
  - `matplotlib`




## Usage

Change parameters in include/parameters.hpp

### Compile

Build, run and plot all problems
```bash
make all
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

### What different scripts do:
- **Example** -

### Plotting
Plot all plots
```bash
make plot
```

