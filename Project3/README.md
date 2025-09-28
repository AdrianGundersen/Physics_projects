# Project 2

## Description


### Project structure
```text
Project3/
├─ include/
│  ├─ constants.hpp
│  ├─ integrator.hpp
│  ├─ parameters.hpp
│  ├─ Particle.hpp
│  └─ PenningTrap.hpp
├─ src/
│  ├─ integrator.cpp
│  ├─ Particle.cpp
│  └─ PenningTrap.cpp
├─ main.cpp
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

Change parameters in include/parameters.hpp

### Compile

Build all problems
```bash
make
```

Clean objects and executables:
```bash
make clean
```
### Run
```bash
From Project root:
./build/bin/project3.exe
```

Or build, then run:
```bash
make run
```

### What different scripts do:
- **Example** -

### Plotting