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
- `problem_2.cpp` - C++ program for Problem 1 & 2 (analytical solution + simple solver)
- `problem_7.cpp` - C++ program implementing the Thomas algorithm (general solver)
- `problem_8.cpp` - C++ program implementing the Thomas algorithm (general solver), finds the error and loos for different iterations $n$
- `problem_10.cpp` - C++ program that comperes the diffrent Thomas algorithms and creates a txt file
- `problem_2_plotting.py` - Python plotting script for Problems 1, 2 & 7 
- `problem_8_plot.py` - Python plotting script for Problem 8
- `output/diff_eq_sol.txt` - Example output file (solution data) are saved in `output` folder
- `filepath.pdf` - Example plot generated from `problem_2_plotting.py`

---

## Usage

### Compile
For Problem `i`:
```bash
g++ -O2 -std=c++17 problem_i.cpp -o problem_i.exe 
```

### Run
For problem 1, 2:
```bash
./problem_i.exe 
```

For problem 7:
```bash
./problem_7.exe filename.txt 
```
- Enter the number of steps $n$ between [0,1] in therminal
- The output file `filename.txt` will  be saved in the `output` folder.
```

For problem 8:
```bash
./problem_8.exe 
```
- Press `N` in therminal. 
- Creates `problem_8_n.txt` files where n in number of iterations and saves in the `output` folder.

For problem 10:
```bash
./problem_10.exe 
```
- Creates `time_optimized_algo.txt` and saves in the `output` folder.

For Python files
```bash
python 3 filename.py
```

### Plot
To plot Problems 1,2 $ 7 run `plotting_sol_1_2.py` 
- Enter the filename in the terminal. It wil autoselect the `output` folder. 
- Select title for the plot. 

To plot problem 8
- Run `plot_probblem_8.py`
