

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


**Example**
How we did the first search for peaks at $L=20$.
```json
{
    "grid" : {
        "M" : 10,
        "L" : 1
    },
    "simulation" : {
        "steps" : 1000,
        "time_step" : 0.01
    },
    "output" : {
        "file_name" : "default"
    },
    "solver" : {
        "method" : "jacobi",
        "tolerance" : 1e-5,
        "max_iterations" : 10000
    },
    "potential" : {
        "type" : "harmonic",
        "frequency" : 1.0
    }
}
```

**Caveats**

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
