# Project 2

## Description

### Project structure
Project2/
├─ include/
│  ├─ jacobi.hpp
│  └─ tridiag.hpp
├─ src/
│  ├─ jacobi.cpp
│  └─ tridiag.cpp
├─ problems/
│  ├─ problem2.cpp
│  ├─ problem3.cpp
│  └─ problem4.cpp
├─ build/                 
├─ Makefile              
└─ README.md

## Usage

### Compile

Build all problems
```bash
make
```
```bash
Or build a specific problem
make problem2      # -> build/problem2.exe
make problem3      # -> build/problem3.exe
make problem4      # -> build/problem4.exe
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
```

### Plot

