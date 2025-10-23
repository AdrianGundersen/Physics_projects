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
```
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
    └── (e.g., analyze.py, plot.ipynb)

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
  - TBD

- **Victor:**  
  - Implement JSON input
  - Grammar, spell correction and proper source documentation
  - Structure folder and skeleton code
  - Make PROJECT_PLAN.md template
  
- **TBD**
  - Optimization of code with OMP and FLOP reduction etc.
  - Report overview and writing
  - Analytical calculations
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
|1 | Make project plan | Everyone | 23.10.2025 | 28.10.2025 | In Progress |
|2 | Make LaTeX template | Victor | 24.10.2025 | 24.10.2025 | In Progress |
| Week 44 |                           |             |            |            |            |
|3 | Problem 1:  2x2 lattice analytics| TBD | TBD | TBD| Not Started |
|4 | Problem 2: Find $\Delta E$ for lattice $L>2$ | TBD | TBD | TBD | Not Started |
|5 | Problem 3: Implement BC without multiple if-tests | TBD | TBD | TBD | Not Started |
|6 | Set up code structure/skeleton code | TBD | 24.10.2025 | TBD | Not Started |
|7 | Implement JSON input template | Victor | 25.10.2025 | 28.10.2025 | Not Started |
---------------------------------------------------------------------------------

## Notes
- Units: 