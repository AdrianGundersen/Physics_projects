# Project 2

## Description

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

### Plot

