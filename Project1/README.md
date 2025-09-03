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
- `output/diff_eq_sol.txt` – Example output file (solution data) are saved in `output` folder
- `Plot2.pdf` – Example plot generated from `plotting_sol_1_2.py`

---

## Usage

### Compile
For Problem 1, 2 & 8:
```bash
g++ -O2 -std=c++17 problem_1_2.cpp -o problem_1_2.exe 
```

For Problem 7:
```bash
g++ -O2 -std=c++17 problem_7.cpp -o problem_7.exe
```

### Run
For problem 1, 2 & 8:
```bash
./problem_1_2.exe 
```

For problem 7:
```bash
./problem_7.exe filename.txt 
```
- Enter the number of steps $n$ between [0,1] in therminal
- The output file `filename.txt` will  be saved in the `output` folder.



### Plot
To plot Problems 1,2 $ 7 run `plotting_sol_1_2.py` 
- Enter the filename in the terminal. It wil autoselect the `output` folder. 
- Select title for the plot. 

To plot problem 8
- Run `problem_8.cpp` first.
- Run `plot_probblem_8.py`
