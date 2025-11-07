# Makefile — FYS3150 Ising project
# Made largely with the help of ChatGPT

SHELL := /bin/bash
.DEFAULT_GOAL := help

# --- Compilers & flags ---
CXX ?= g++
MODE ?= release          # release | debug
ARCH ?=                  # e.g. -march=native or Pi flags
CXXFLAGS := -std=c++20 -fopenmp $(ARCH)
ifeq ($(MODE),release)
CXXFLAGS += -O3
else ifeq ($(MODE),debug)
CXXFLAGS += -O0 -g -DDEBUG
endif

INCLUDE := -Iinclude -I$(abspath ..)/third_party

# --- Paths ---
BIN_DIR := bin

# --- Core sources ---
SRC_CORE := \
  src/ising/lattice.cpp \
  src/ising/observables.cpp \
  src/ising/metropolis.cpp \
  src/ising/mcmc_run.cpp \
  src/omp_rng.cpp

# --- Apps ---
APP_VALIDATE := apps/validate2x2.cpp
BIN_VALIDATE := $(BIN_DIR)/validate2x2
JSON_2x2     ?= configs/2x2_test.json

APP_L20 := apps/L20.cpp
BIN_L20 := $(BIN_DIR)/L20
JSON_L20 ?= configs/L20.json

APP_LN := apps/Ln.cpp
BIN_LN := $(BIN_DIR)/Ln
JSON_LN ?= configs/multiple_walkers.json

APP_RUNTIME := apps/runtime.cpp
BIN_RUNTIME := $(BIN_DIR)/runtime
JSON_RUNTIME ?= configs/multiple_walkers.json

# --- Aggregate ---
APPS := $(BIN_VALIDATE) $(BIN_L20) $(BIN_LN) $(BIN_RUNTIME)

# --- Build rules (single-step link, no object dir) ---
$(BIN_DIR):
	@mkdir -p $@

$(BIN_VALIDATE): $(SRC_CORE) $(APP_VALIDATE) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@

$(BIN_L20): $(SRC_CORE) $(APP_L20) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@

$(BIN_LN): $(SRC_CORE) $(APP_LN) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@

$(BIN_RUNTIME): $(SRC_CORE) $(APP_RUNTIME) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@

# --- Top-level ---
all: $(APPS)

run-2x2: $(BIN_VALIDATE)
	./$(BIN_VALIDATE) $(JSON_2x2)

run-L20: $(BIN_L20)
	./$(BIN_L20) $(JSON_L20)

run-LN: $(BIN_LN)
	./$(BIN_LN) $(JSON_LN)

run-runtime: $(BIN_RUNTIME)
	./$(BIN_RUNTIME) $(JSON_RUNTIME)

run: run-LN

# --- Utilities ---
format:
	@cf=$$(command -v clang-format || true); \
	test -n "$$cf" || { echo "clang-format not found"; exit 1; }; \
	find include src apps -type f \( -name '*.hpp' -o -name '*.h' -o -name '*.cpp' \) -print0 | xargs -0 $$cf -i

clean:
	rm -rf $(BIN_DIR)

help:
	@echo "Targets:"
	@echo "  all                Build all apps (MODE=$(MODE))"
	@echo "  run-2x2            Run validate2x2 (JSON_2x2)"
	@echo "  run-L20            Run L20 (JSON_L20)"
	@echo "  run-LN             Run Ln  (JSON_LN)"
	@echo "  run-runtime        Run runtime (JSON_RUNTIME)"
	@echo "  format             clang-format sources"
	@echo "  clean              Remove bin/"
	@echo
	@echo "Remote Raspberry Pi:"
	@echo "  pirsync            rsync project to $$PI_HOST:$$PI_DIR (skips bin/.git)"
	@echo "  pibuild            build remotely with ARCH='$$PI_ARCH' MODE='$$MODE'"
	@echo "  run-pi-2x2         build on Pi, then run validate2x2"
	@echo "  run-pi-L20         build on Pi, then run L20"
	@echo "  run-pi-LN          build on Pi, then run Ln"
	@echo
	@echo "Variables:"
	@echo "  MODE=release|debug       (default: release)"
	@echo "  ARCH='-march=native'     (empty by default; Pi uses PI_ARCH)"
	@echo "  JSON_2x2 / JSON_L20 / JSON_LN / JSON_RUNTIME"

.PHONY: all run-2x2 run-L20 run-LN run-runtime format clean help

# --- Pi targets ---
PI_HOST ?= rpi
PI_DIR  ?= ~/$(notdir $(CURDIR))
PI_ARCH ?= -march=armv8-a+simd -mcpu=cortex-a76
RSYNC_EXCLUDES := --exclude bin --exclude .git --exclude __pycache__

pirsync:
	rsync -az $(RSYNC_EXCLUDES) ./ $(PI_HOST):$(PI_DIR)/

pibuild: pirsync
	ssh $(PI_HOST) 'cd $(PI_DIR) && $(MAKE) MODE=$(MODE) ARCH="$(PI_ARCH)" -j$$(nproc)'

run-pi-2x2: pibuild
	ssh $(PI_HOST) 'cd $(PI_DIR) && ./$(BIN_VALIDATE) $(JSON_2x2)'

run-pi-L20: pibuild
	ssh $(PI_HOST) 'cd $(PI_DIR) && ./$(BIN_L20) $(JSON_L20)'

run-pi-LN: pibuild
	ssh $(PI_HOST) 'cd $(PI_DIR) && ./$(BIN_LN) $(JSON_LN)'

.PHONY: pirsync pibuild run-pi-2x2 run-pi-L20 run-pi-LN