# Project 2

## Description
This project implements and analyzes the Jacobi method for computing eigenvalues and eigenvectors of real symmetric matrices, with the buckling-beam problem as the primary application. We compare numerical results against analytical expressions for the tridiagonal stiffness matrix that arises from a second-order finite-difference discretization.

### Project structure
```text
Project2/
├─ include/
│  ├─ jacobi.hpp
│  ├─ tridiag.hpp
│  └─ analytical.hpp
├─ src/
│  ├─ jacobi.cpp
│  ├─ tridiag.cpp
│  └─ analytical.cpp
├─ problems/
│  ├─ problem2.cpp
│  ├─ problem3.cpp
│  ├─ problem4.cpp
│  ├─ problem5.cpp
│  └─ problem6.cpp
├─ build/
├─ Makefile
└─ README.md
```

## Dependencies

- **C++20 compiler** (should also work down to C++17)  
- **Armadillo**   
- **Python 3.10+**
  - `numpy`
  - `matplotlib`



## Usage

### Compile

Build all problems
```bash
make
```
Or build a specific problem
```bash
make problem2      # -> build/problem2.exe
make problem3      # -> build/problem3.exe
make problem4      # -> build/problem4.exe
make problem5      # -> build/problem5.exe
make problem6      # -> build/problem6.exe

```

Clean objects and executables:
```bash
make clean
```
### Run
```bash
From Project root:
./build/problem2.exe
./build/problem3.exe
./build/problem4.exe
./build/problem5.exe
./build/problem6.exe
```

Or build, then run:
```bash
make run2
make run3
make run4
make run5
make run6
```

### What each problem does:
- **Problem 2** - Tests armadillos eigenvalue and eigenvector functions against known analytical values. 
- **Problem 3** - Finds maximum off-diagonal matrix-element of a symmetrical matrix.
- **Problem 4** - Implement an algorithm to compute eigenvalues and eigenvactors of a $N\times N$ tridiagonal matrix with $N=6$ satisfying $\textbf{A} \vec{v} = \lambda \vec{v}$, and compare the result to the analytical solution.
- **Problem 5** - Runs the algorithm in problem 4 for different values of N to estimate how the number of iterations scale. When runs the algorithm for a dense symetric matrix to compare the scaling. 
- **Problem 6** - Runs the algorithm for $N=9$ and $N=99$. Adds boundery points $(\hat{x}_0, \hat{v}_0)$ and $(\hat{x}_n, \hat{v}_n)$. Plots the eigenvectors whom corresponds to the three lowest eigenvalues against $\hat{x}\in [0,1]$. Repeats for the analitical solution for comparison. 

### Plotting
Problem 5:
```bash
python3 /problems/problem5.py
```

Problem 6:
```bash
python3 /problems/problem6.py
```
