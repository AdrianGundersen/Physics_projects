# Makefile (project root)
CXX := g++
ARCH ?=
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra $(ARCH)

INCLUDE := -Iinclude -I$(abspath ..)/third_party # third party at FYS3150 root
BIN := bin/test_json 
SRC := test/JSON_example/testJSON.cpp
# App: validate 2x2
BIN_DIR := bin
BIN_VALIDATE := $(BIN_DIR)/validate2x2
APP_VALIDATE := apps/validate2x2.cpp
SRC_CORE := src/ising/lattice.cpp src/ising/observables.cpp src/ising/metropolis.cpp    # add more .cpp later as you implement

all: $(BIN_VALIDATE)

$(BIN_VALIDATE): $(SRC_CORE) $(APP_VALIDATE)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@

# Default config to run (2x2)
JSON ?= configs/2x2_test.json

run-2x2: $(BIN_VALIDATE)
	./$(BIN_VALIDATE) $(JSON) 

clean:
	rm -rf $(BIN_DIR)

.PHONY: all run-2x2 clean

# --- Pi remote build (optional) ---
PI_HOST ?= rpi
PI_DIR  ?= ~/$(notdir $(CURDIR))
PI_ARCH ?= -march=armv8-a+simd -mcpu=cortex-a76

pirsync:
	rsync -az --delete ./ $(PI_HOST):$(PI_DIR)/

pibuild: pirsync
	ssh $(PI_HOST) 'cd $(PI_DIR) && \
		$(MAKE) clean || true; \
		$(MAKE) ARCH="$(PI_ARCH)" -j$$(nproc)'

run-pi: pibuild
	ssh $(PI_HOST) 'cd $(PI_DIR) && ./$(BIN_VALIDATE) $(JSON)'

.PHONY: pirsync pibuild run-pi
