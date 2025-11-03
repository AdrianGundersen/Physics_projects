# Project Plan
A description of the project, goals, timeline, and task assignments. Only meant to be used internally.

## Project Overview

**Ising model simulation** on an $L \times L$ lattice using the Metropolis algorithm.

### Milestones and Primary Goals
- **Analytical calculations**
- **Core code** structure and JSON input
- **Optimization** with OpenMP
- **Analysis and plotting** 
- **Report writing** and review
- **Final delivery** of code and report
### File Structure
```text
Project4/
├── Makefile                 # Build targets (consider switching to CMake if desired)
├── README.md                # Project overview and quickstart
├── PROJECT_PLAN.md          # Internal project plan (scope, timeline, tasks)
│
├── configs/                 # JSON configs for different runs
│   └── 2x2_test.json
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
│   │   ├── io/
│   │   │   └── json_util.hpp # JSON parsing utilities (Header-only)
│   │   ├── lattice.hpp       # Lattice class definition with L and spins
│   │   └── metropolis.hpp    # Metropolis algorithm functions
│   │   └── observables.hpp   # Functions to compute observables like energy, magnetization
│   │   └── model.hpp         # Model parameters and model structure (Header-only)
│   └── omp_rng.hpp       # OpenMP random number generator
│
└── src/
    └─ ising/
    |  ├─ lattice.cpp
    |  ├─ metropolis.cpp
    |  └─ observables.cpp
    └─ omp_rng.cpp
│
├── apps/                    # Executables running JSON configs
│   └── validate2x2.cpp
│
└── scripts/                 # Analysis & plotting (Python scripts)
    └── 
```
## Timeline
Project Start Date: 23.10.2025
Project End Date: 19.11.2025

**Week 43** *(Thursday 23.10.2025 - Sunday 26.10.2025)*
- Make draft for project plan

**Week 44** *(Monday 27.10.2025 - Sunday 02.11.2025)*
- Do all analytical calculations
- Set up basic code structure/skeleton code
- Implement JSON input template for easy runs
  
**Week 45** *(Monday 03.11.2025 - Sunday 09.11.2025)*
- Parallelization with OMP
  
**Week 46** *(Monday 10.11.2025 - Sunday 16.11.2025)*
- Be done with all analysis and plotting
- Write report
  
**Week 47** *(Monday 17.11.2025 - Tuesday 18.11.2025)*
- Grammar and spell check of report
- Review code and documentation
- Deliver 19.11.2025
  
## Tasks and Assignments

### Main responsibilities
*Note*: This is just a preliminary distribution of tasks. Everyone is expected to keep track of all parts.

Things to pick:

- **Adrian:**
  - GitHub management
  
- **Casper:** 
  - Analytical calculations for 2x2 lattice

- **Victor:**  
  - Implement JSON input
  - Grammar, spell correction and proper source documentation
  - Structure folder and skeleton code
  - Make PROJECT_PLAN.md template
  
- **TBD**
  - Optimization of code with OMP and FLOP reduction etc.
  - Report overview and writing

  - Plotting and analysis
  - Makefile
  - Readme
  - LaTeX template
  - JSON input implementation
  - Code documentation
---------------------------------------------------------------------------------

### Task List
**Note**: Should weekly agree on tasks to be done and update the plan according to status.
| Task ID | Task Description      | Assigned To | Start Date | End Date   | Status     |
|---------|---------------------------|-------------|------------|------------|------------|
| Week 43 |                           |             |            |            |            |
|1 | Make project plan draft | Everyone | 23.10.2025 | 28.10.2025 | Completed |
|2 | Make LaTeX template | Victor | 24.10.2025 | 24.10.2025 | Completed |
| Week 44 |                           |             |            |            |            |
|3 | Problem 1:  2x2 lattice analytical | Casper & Victor | 25.10.2025 | 30.10.2025 | Completed |
|4 | Problem 2: Find $\Delta E$ for lattice $L>2$ | Victor | 26.10.2025 | 26.10.2025 | Completed |
|5 | Problem 3: Implement BC without multiple if-tests | Adrian | 27.10.2025 | 28.10.2025 | Completed |
|6 | Set up code structure/skeleton code | Victor | 24.10.2025 | 30.10.2025 | Completed |
|7 | Implement JSON input template | Victor | 25.10.2025 | 28.10.2025 | Completed |
| 8| Add option to run on raspberry pi | Victor | 26.10.2025 | 26.10.2025 | Completed |
|9 | Calculate observables| Adrian & Victor | 29.10.2025 | 29.10.2025 | Completed |
|10 | Write Boltzmann lookup table | Victor | 29.10.2025 | 30.10.2025 | Completed |
|11 | Implement Metropolis algorithm | Adrian & Casper | 30.10.2025 | 02.11.2025 | Completed |
|12 | Validate 2x2 lattice code | Adrian & Casper | 01.11.2025 | 02.11.2025 | Completed |
|13 | Study burn-in time| Adrian | 01.11.2025 | 02.11.2025 | Completed |
|14 | Find probability distribution of states | Adrian & Victor | 02.11.2025 | 02.11.2025 | Completed |
| Week 45 |                           |             |            |            |            |
|15 | Parallelize code with OpenMP as a MCMC_run.cpp | Adrian | 03.11.2025 | 03.11.2025 | Completed |
|16 | Calculate speedup | Victor | 03.11.2025 | 03.11.2025 | Completed |
| Week 46 |                           |

## Notes
- Units: 