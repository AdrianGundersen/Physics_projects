# Project 2

## Description


## Usage

### Compile

Shared object files:
```bash
g++ -O2 -std=c++20 -c tridiag.cpp -o tridiag.o
g++ -O2 -std=c++20 -c jacobi.cpp -o jacobi.o
```
For Problem `2`:
```bash
g++ -O2 -std=c++20 problem2.cpp tridiag.o -larmadillo -o problem2.exe
```
For Problem `3`:
```bash
g++ -O2 -std=c++20 problem3.cpp tridiag.o -larmadillo -o problem3.exe
```

For problem '4'
```bash
g++ problem4.cpp jacobi.cpp tridiag.cpp -o problem4.exe -larmadillo
```

### Run


### Plot

