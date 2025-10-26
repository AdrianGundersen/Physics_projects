# Makefile (project root)
CXX := g++
# Optional architecture flags (to use raspberry pi)
ARCH ?=

CXXFLAGS := -std=c++20 -O2 -Wall -Wextra $(ARCH)

INCLUDE := -Iinclude -I$(abspath ..)/third_party

BIN := bin/test_json
SRC := test/JSON_example/testJSON.cpp

all: $(BIN)

$(BIN): $(SRC)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -o $@

JSON ?= test/JSON_example/testJSON.json

run: $(BIN)
	./$(BIN) $(JSON)

clean:
	rm -rf bin

.PHONY: all run clean


# Pi remote build (optional)
# SSH alias or user@host of your Pi
PI_HOST ?= rpi
# Remote directory to sync into (uses the current folder name)
PI_DIR  ?= ~/$(notdir $(CURDIR))

# CPU tuning (Pi 5 default)
PI_ARCH ?= -march=armv8-a+simd -mcpu=cortex-a76

.PHONY: pirsync pibuild run-pi

pirsync:
	rsync -az --delete ./ $(PI_HOST):$(PI_DIR)/

pibuild: pirsync
	ssh $(PI_HOST) 'cd $(PI_DIR) && \
		$(MAKE) clean || true; \
		$(MAKE) ARCH="$(PI_ARCH)" -j$$(nproc)'

run-pi: pibuild
	ssh $(PI_HOST) 'cd $(PI_DIR) && ./$(BIN) $(JSON)'
