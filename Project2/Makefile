# Made largely by ChatGPT, with slight modifications
# ---------- Compiler & flags ----------
CXX      := g++
CXXFLAGS := -O2 -std=c++20 -Iinclude
LDFLAGS  := -larmadillo

# ---------- Folders ----------
SRC_DIR   := src
PROB_DIR  := problems
BUILD_DIR := build

# ---------- Library sources/objects ----------
LIB_SRCS  := $(SRC_DIR)/tridiag.cpp $(SRC_DIR)/jacobi.cpp
LIB_OBJS  := $(addprefix $(BUILD_DIR)/,$(notdir $(LIB_SRCS:.cpp=.o)))

# ---------- Problem sources ----------
P2_SRCS := $(PROB_DIR)/problem2.cpp
P3_SRCS := $(PROB_DIR)/problem3.cpp
P4_SRCS := $(PROB_DIR)/problem4.cpp
P5_SRCS := $(PROB_DIR)/problem5.cpp
P6_SRCS := $(PROB_DIR)/problem6.cpp

# ---------- Executables ----------
P2_EXE := $(BUILD_DIR)/problem2.exe
P3_EXE := $(BUILD_DIR)/problem3.exe
P4_EXE := $(BUILD_DIR)/problem4.exe
P5_EXE := $(BUILD_DIR)/problem5.exe
P6_EXE := $(BUILD_DIR)/problem6.exe

.PHONY: all dirs clean problem2 problem3 problem4 problem5 problem6 run2 run3 run4 run5 run6 help

all: dirs $(P2_EXE) $(P3_EXE) $(P4_EXE) $(P5_EXE) $(P6_EXE)

dirs:
	@mkdir -p $(BUILD_DIR)

# Build library objects from src/
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | dirs
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---------- Link rules per problem ----------
# 2 & 3 use tridiag only
$(P2_EXE): $(P2_SRCS) $(BUILD_DIR)/tridiag.o | dirs
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(P3_EXE): $(P3_SRCS) $(BUILD_DIR)/tridiag.o | dirs
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# 4, 5, 6 use jacobi + tridiag
$(P4_EXE): $(P4_SRCS) $(BUILD_DIR)/jacobi.o $(BUILD_DIR)/tridiag.o | dirs
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(P5_EXE): $(P5_SRCS) $(BUILD_DIR)/jacobi.o $(BUILD_DIR)/tridiag.o | dirs
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(P6_EXE): $(P6_SRCS) $(BUILD_DIR)/jacobi.o $(BUILD_DIR)/tridiag.o | dirs
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# ---------- Convenience aliases ----------
problem2: $(P2_EXE)
problem3: $(P3_EXE)
problem4: $(P4_EXE)
problem5: $(P5_EXE)
problem6: $(P6_EXE)

run2: $(P2_EXE)
	$<

run3: $(P3_EXE)
	$<

run4: $(P4_EXE)
	$<

run5: $(P5_EXE)
	$<

run6: $(P6_EXE)
	$<

clean:
	@rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.exe

help:
	@echo "Targets:"
	@echo "  make              -> build all problems (2..6)"
	@echo "  make problemN     -> build only problemN.exe (N=2..6)"
	@echo "  make runN         -> build & run problemN.exe (N=2..6)"
	@echo "  make clean        -> remove objects and executables"
