# Makefile — FYS3150 Project 5 Schrödinger CN solver

SHELL := /bin/bash
.DEFAULT_GOAL := help

# --- Compiler & flags ---
CXX  ?= g++
MODE ?= release          # release | debug
ARCH ?=                  # e.g. -march=native

CXXFLAGS := -std=c++20 -Wall -Wextra $(ARCH)
ifeq ($(MODE),release)
CXXFLAGS += -O3
else ifeq ($(MODE),debug)
CXXFLAGS += -O0 -g -DDEBUG
endif

INCLUDE := -Iinclude -I$(abspath ..)/third_party
LDFLAGS :=

# --- Paths ---
BIN_DIR := bin
OUT_DIR := data/outputs

# --- Sources ---
SRC_CORE := \
  src/solver.cpp \
  src/io.cpp \
  src/potential.cpp

APP_MAIN := src/main.cpp
BIN_MAIN := $(BIN_DIR)/ds_solver
CONFIG   ?= configs/configs.json

# --- Aggregate ---
APPS := $(BIN_MAIN)

# --- Build rules (single-step link, no object dir) ---
$(BIN_DIR):
	@mkdir -p $@

$(OUT_DIR):
	@mkdir -p $@

$(BIN_MAIN): $(SRC_CORE) $(APP_MAIN) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LDFLAGS)

# --- Top-level ---
all: $(APPS) | $(OUT_DIR)

run: $(BIN_MAIN)
	./$(BIN_MAIN) $(CONFIG)

# --- Utilities ---
clean:
	rm -rf $(BIN_DIR)

help:
	@echo "Targets:"
	@echo "  all           Build ds_solver (MODE=$(MODE))"
	@echo "  run           Run ds_solver with CONFIG"
	@echo "  clean         Remove bin/"
	@echo
	@echo "Variables:"
	@echo "  MODE=release|debug          (default: release)"
	@echo "  ARCH='-march=native'        (empty by default)"
	@echo "  CONFIG=path/to/config.json  (default: configs/configs.json)"

.PHONY: all run clean help
