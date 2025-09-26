# ---------- Compiler & flags ----------
CXX      := g++
CXXFLAGS := -O2 -std=c++20 -Iinclude
LDFLAGS  := -larmadillo

# ---------- Folders ----------
SRC_DIR   := src
BUILD_DIR := build

# ---------- Sources & objects ----------
SRCS := $(SRC_DIR)/Particle.cpp $(SRC_DIR)/PenningTrap.cpp $(SRC_DIR)/main.cpp
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# ---------- Executable ----------
EXE := $(BUILD_DIR)/project3.exe

.PHONY: all clean run dirs

all: dirs $(EXE)

dirs:
	@mkdir -p $(BUILD_DIR)

# Compile object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | dirs
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link final executable
$(EXE): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Convenience run target
run: $(EXE)
	$<

clean:
	@rm -rf $(BUILD_DIR)
