# Made largely by ChatGPT, with slight modifications

# ---------- Compiler & flags ----------
CXX      := g++
CXXFLAGS := -O2 -std=c++20 -Iinclude
LDFLAGS  := -larmadillo

# ---------- Folders ----------
SRC_DIR   := src
PROB_DIR  := problems
BUILD_DIR := build

# ---------- Library objects ----------
LIB_SRCS  := $(SRC_DIR)/tridiag.cpp $(SRC_DIR)/jacobi.cpp
LIB_OBJS  := $(addprefix $(BUILD_DIR)/,$(notdir $(LIB_SRCS:.cpp=.o)))

# ---------- Problem sources ----------
P2_SRCS := $(PROB_DIR)/problem2.cpp
P3_SRCS := $(PROB_DIR)/problem3.cpp
P4_SRCS := $(PROB_DIR)/problem4.cpp

# ---------- Executables ----------
P2_EXE := $(BUILD_DIR)/problem2.exe
P3_EXE := $(BUILD_DIR)/problem3.exe
P4_EXE := $(BUILD_DIR)/problem4.exe

.PHONY: all dirs clean run2 run3 run4 help

# Build everything by default
all: dirs $(P2_EXE) $(P3_EXE) $(P4_EXE)

# Ensure build dir exists
dirs:
	@mkdir -p $(BUILD_DIR)

# Pattern rule for library object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | dirs
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---------- Link rules per problem ----------
# Problem 2 uses tridiag only
$(P2_EXE): $(P2_SRCS) $(BUILD_DIR)/tridiag.o | dirs
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Problem 3 uses tridiag only
$(P3_EXE): $(P3_SRCS) $(BUILD_DIR)/tridiag.o | dirs
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Problem 4 uses jacobi + tridiag
$(P4_EXE): $(P4_SRCS) $(BUILD_DIR)/jacobi.o $(BUILD_DIR)/tridiag.o | dirs
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# ---------- Convenience run targets ----------
run2: $(P2_EXE)
	$<

run3: $(P3_EXE)
	$<

run4: $(P4_EXE)
	$<

# ---------- Cleanup ----------
clean:
	@rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.exe

# ---------- Help ----------
help:
	@echo "Targets:"
	@echo "  make            -> build all problems (2,3,4)"
	@echo "  make problem2   -> build only problem2.exe"
	@echo "  make problem3   -> build only problem3.exe"
	@echo "  make problem4   -> build only problem4.exe"
	@echo "  make run2/run3/run4 -> build & run selected problem"
	@echo "  make clean      -> remove objects and executables"
