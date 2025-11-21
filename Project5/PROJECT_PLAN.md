# Project Plan
A description of the project, goals, timeline, and task assignments. Only meant to be used internally.

## Project Overview
**Double slit experiment simulation** on a 2D grid using the Crank-Nicolson method to solve the time-dependent Schrödinger equation.

Linear system:

$$
A \mathbf{v}^{n+1} = B \mathbf{v}^{n}
$$

however we never form the matrices explicitly, but use stencil formulas to access neighboring points in the grid.


## Useful classes, algorithms, and methods to implement

- We will not be using matrices but rather complex arrays from the C++ standard library. NO ARMADILLO!!! Slow piece of junk - I need big grid for pretty solution!
- Instead of storing big matrices we use stencil formulas to access neighboring points in the grid (5 point stencil).

- constants
   - physical constants = 1
   - numerical constants = grid size, time step, total time, number of time steps

- arrays
    - complex RHS
    - complex precomputed coefficients for inverse matrix
    - v_old wavefunction (size M*M) (assume square grid)
    - v_new wavefunction (size M*M)
  

- build_rhs 
   - constructs the right hand side 
   - input: v_old
   - output: rhs
   - uses stencil formula to access neighboring points
   - applies potential

- precompute_coeff
    - precomputes coefficients for the matrix
    - Dij (diagonal elements of A matrix)
    - inv_Dij (inverse of diagonal elements)
    - used in Jacobi algorithm

- JacobiSolver
    - solves the system using the Jacobi algorithm
    - guess initial solution
    - find new solution iteratively until convergence

- Crank-Nicolson step
  - performs a single time step using the Crank-Nicolson method
  - uses Jacobi algorithm to solve the system
  - precomputed coefficients for the matrix

### Milestones and Primary Goals

### File Structure
```text
Project5/
├── Makefile                     # Build targets for the C++ code
├── README.md                    # Project overview and quickstart
├── PROJECT_PLAN.md              # Internal project plan 
├── include/                     # C++ headers
│   ├── constants.hpp            # Physical and numerical constants
│   ├── io.hpp                   # Input/output function declarations
│   ├── solver.hpp               # Crank-Nicolson + Jacobi
│   └── potential.hpp            # Functions to build the double-slit potential and walls
├── src/                         # C++ source files
│   ├── main.cpp                 # Sets up simulation, parses config, runs time loop
│   ├── solver.cpp               # Implementation of solver routines
│   ├── potential.cpp            # Implementation of potential construction
│   └── io.cpp                   # Implementation of I/O utilities
|
├── configs/                     # Config files for different runs
│   └── ...
|
├── data/                        # Outputs
│   └── ... 
|
├── scripts/                     # Analysis & plotting scripts (Python)
│   ├── plot.py                  # Plot single snapshots
│   └── animate.py               # Create animations of the double-slit
|
└── figs/                        # Generated figures and animations
    └── 


```
## Timeline
Project Start Date: 21.11.2025
Project End Date: 10.12.2025

**Week 47** *(Friday 21.11.2025 - Sunday 23.11.2025)*
- 
**Week 48** *(Monday 24.11.2025 - Sunday 30.11.2025)*

**Week 49** *(Monday 01.12.2025 - Sunday 07.12.2025)*

**Week 50** *(Monday 08.12.2025 - Wednesday 10.12.2025)*

  
## Tasks and Assignments

### Main responsibilities
*Note*: This is just a preliminary distribution of tasks. Everyone is expected to keep track of all parts.

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
---------------------------------------------------------------------------------

### Simulation Plan



### Task List
**Note**: Should weekly agree on tasks to be done and update the plan according to status.
| Task ID | Task Description      | Assigned To | Start Date | End Date   | Status     |
|---------|---------------------------|-------------|------------|------------|------------|
| Week 47 |                           |             |            |            |            |
|1 | Make project plan draft | Everyone | 23.10.2025 | 28.10.2025 | In Progress |   
| 2 | Make LaTeX template | Victor | 24.10.2025 | 24.10.2025 | Completed |

## Notes
- Units: 