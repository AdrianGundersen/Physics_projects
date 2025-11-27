

# Project 5 - The Time-Dependent Schrödinger Equation in 2D & The Double Slit Experiment

## Description
We simulate the time-dependent Schrödinger equation on a 2D grid using the Crank-Nicolson method. The simulation focuses on the double slit experiment, demonstrating wavefunction evolution and interference patterns. 

One of our main focuses is on how to implement Crank-Nicolson without forming large matrices and doing heavy linear algebra, but use stencil formulas. 


## Physical Model
The time-dependent Schrödinger equation in two dimensions is given by:

$$
\mathrm{i}\hbar \frac{\partial}{\partial t} \ket{\Psi}(x,y,t) = H\ket{\Psi}
$$

$$
\mathrm{i}\hbar \frac{\partial \Psi(x,y,t)}{\partial t} = -\frac{\hbar^2}{2m} \left( \frac{\partial^2 \Psi(x,y,t)}{\partial x^2} + \frac{\partial^2 \Psi(x,y,t)}{\partial y^2} \right) + V(x,y) \Psi(x,y,t)
$$

or in dimensionless units:

$$
u_t = \mathcal{L}u
$$

where

 $$
\mathcal{L} = \mathrm{i}\Delta - \mathrm{i} V(x,y)
$$


## Project Structure
```text

```

## Plot Config overview


## Installation and Dependencies


- **C++20 compiler** (e.g. g++ 11 or higher)
- **Armadillo** (Not used for the main solver, but for solving matrix problems as requested per the project description)
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


**Example**
How we did the first search for peaks at $L=20$.
```json
{
    "grid" : {
        "M" : 200,          # number of grid points in x and y
        "L" : 1             # length of box domain
    },
    "simulation" : {
        "N" : 81,           # total number of time steps t0, ... tN
        "dt" : 2.5e-5,      # time step size
        "xC" : 0.25,        # x center of initial wave packet
        "sigma_x" : 0.05,   # standard deviation in x
        "px" : 200,         # initial momentum in x (wavenumber for unitless)
        "yC" : 0.5,         # y center of initial wave packet
        "py" : 0,           # initial momentum in y (wavenumber for unitless)
        "sigma_y" : 0.2     # standard deviation in y

    },
    "output" : {
        "file_name" : "output/probability_density_p8.txt",                  # file for probability density output (not used, to use uncomment in solver.cpp)
        "file_name_wavefunction" : "output/wavefunction_single_slit.txt",   # file for wavefunction output    
        "precision" : 8                                                     # output precision
    },
    "solver" : {
        "method" : "jacobi",     # solver method (only jacobi/gauss_seidel implemented)
        "tolerance" : 1e-4,      # solver tolerance for jacobi
        "max_iterations" : 2000  # max iterations for jacobi
    },
    "potential" : {
        "type" : "none",         # type of base potential (none, harmonic)
        "frequency" : 1.0,       # frequency for harmonic potential
        "V0" : 1e10              # potential strength for double slit wall
    },
    "slits" : { 
        "enabled" : true,        # enable/disable slits
        "wall_center" : 0.5,     # center position of wall in x
        "wall_thickness" : 0.02, # thickness of wall in x
        "slit_aperture" : 0.05,  # aperture of slits
        "num_slits" : 1,         # number of slits
        "slit_spacing" : 0.05    # spacing between slits
    },
    "parallelization" : {
        "threads" : 12          # number of OpenMP threads
    }
}
```

**Caveats**

**Simulation**
- Describes initial gaussian wavepacket. 
- Only Dirichlet BC implemented.

**Potential**
- Harmonic potential was only used for testing, so may not be fully functional.

**Solver**
- Only the Jacobi iterative method is implemented. 
- Slow convergence for large grids. 
- Gives larger errors at first time steps. 
- max_iters is per time step.
- tolerance is the maximum change at any grid point.

**Slits**
- The slits are spaced uniformly around the center of the wall with equal spacing.
- The total width of all slits and spacings must fit within the wall thickness. There is no test for this, so be careful.
- If the wavepacket is initialized over the wall, unexpected behavior may occur.
  

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


### Plotting

### What different scripts do

**Simulation apps**

**Plotting and analysis utilities (Python)**


## Important parameters for reproducability


## Future work

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
