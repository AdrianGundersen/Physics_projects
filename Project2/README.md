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
- **Problem 2** -- Tests armadillos eigenvalue and eigenvector functions against known analytical values. 
- **Problem 3** -- Finds maximum off-diagonal matrix-element of a symmetrical matrix.
- **Problem 4** -- Creates and runs an algorythm to find eigenvalues and eigenvactors of a $N\times N$ tridiagonal matrix where $N=6$ satisfying $\textbf{A} \vec{v} = \lambda \vec{v}$ and compares it to the analytical result.
- **Problem 5** -- Runs the algorythm in problem 4 for different values of N to estimate how the amound of iterations scale. When runs the algorythm for a dense symetric matrix to compare its scaling. 

### Plotting
Problem 5:
```bash
python3 /problems/problem5.py
```

Problem 6:
```bash
python3 /problems/problem6.py
```
