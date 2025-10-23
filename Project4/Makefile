# Makefile (project root)
CXX := g++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra
INCLUDE := -Iinclude

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
