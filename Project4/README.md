# Project 4

## Description

## Physical Model


---


## Project Structure
```text
Project4/
├── Makefile                 # Build targets (consider switching to CMake if desired)
├── README.md                # Project overview and quickstart
├── PROJECT_PLAN.md          # Internal project plan (scope, timeline, tasks)
│
├── configs/                 # JSON configs for different runs
│   └── ...
│
│
├── data/                 # Numerical outputs, logs, tables
│   └── auto-generated data files
│
├── figs/                    # Plots and figures generated from results
│   └── (auto-generated images)
│
├── include/                 # Header files
│   ├── ising/
│   │   ├── lattice.hpp
│   │   └── metropolis.hpp
│   └── constants.hpp
│
├── src/                     # Implementation files
│   └── metropolis.cpp
│
├── apps/                    # Executables running JSON configs
│   └── 
│
└── scripts/                 # Analysis & plotting (Python scripts)
    └── 
```

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
The executables in the `apps/` folder read parameters from JSON config files located in the `configs/` folder. You can create your own config files by copying and modifying the existing templates.

**Example**
```json
{
    "model": // Model parameters
    {
        "J": 1.0, // Interaction strength
        "double_count" : false, // Whether to double count energy
        "spin_config": "all_up" // Initial spin configuration: "all_up", "all_down", "random"
    },
    "lattice": // Lattice parameters
    {
        "L" : 2 // Lattice size (LxL)
    },
    "simulation": // Simulation parameters (Metropolis)
    {
        "seed": 67, // Mother seed for random number generator
        "temperature": 2.0, // Temperature [J/kB]
        "total_steps": 1000, // Total Monte Carlo steps
        "measure_step": 10 // Measure observables every n steps
    }
}
```

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
make run-2x2 JSON=configs/2x2_test.json
```


On raspberry pi, use:
```bash
make PI_HOST=<your_pi_host>
```

or set the `PI_HOST` variable in the `Makefile` and run:
```bash
make run-pi
```


### Plotting

### Example usage:


### What different scripts do:


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