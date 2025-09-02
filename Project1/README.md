# Project 1: 1D Poisson Solver

## Description
Numerical solution of the boundary value problem:

$$
\begin{cases}
    -u''(x) = f(x), & x \in [0,1], \\
    f(x) = 100 e^{-10x}, & \\
    u(0) = 0, \quad u(1) = 0. &
\end{cases}
$$

- Derive the analytical solution
- Discretize with finite differences
- Solve the tridiagonal system using Thomas algorithm
- Compare analytical vs. numerical solutions
- Compute and tabulate maximum relative errors

## Files
- `problem_1_2.cpp` – C++ program for Problem 1 & 2 (analytical solution + simple solver)
- `problem_7.cpp` – C++ program implementing the Thomas algorithm (general solver)
- `plotting_sol_1_2.py` – Python plotting script for Problems 1, 2 & 7 
- `diff_eq_sol.txt` – Example output file (solution data)
- `Plot2.pdf` – Example plot generated from `plotting_sol_1_2.py`

---

## Usage

### Compile
For Problem 1 & 2:
```bash
g++ -O2 -std=c++17 problem_1_2.cpp -o problem_1_2.exe 
```

For Problem 7:
```bash
g++ -O2 -std=c++17 problem_7.cpp -o problem_7.exe
```

### Run
For problem 7:
```bash
./problem_7.exe filename.txt 
```
Enter the number of steps between [0,1] 

### Plot
To plot Problems 1,2 $ 7 run `plotting_sol_1_2.py` 
- Enter the filename in the terminal. It wil autoselect the `output` folder. 
- Select title for the plot. 

