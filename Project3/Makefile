# ---------- Compiler & flags ----------
CXX       ?= g++
CXXSTD     = -std=c++20
INC        = -Iinclude
LIBS       = -larmadillo
# If your toolchain needs it for std::filesystem, uncomment:
# LIBS      += -lstdc++fs

# Build type (make MODE=debug | release | sanitize)
MODE      ?= release
COMMON    = $(CXXSTD) $(INC) -MMD -MP
CXXFLAGS_debug    = -O0 -g3 -fno-omit-frame-pointer
CXXFLAGS_release  = -O3 -march=native
CXXFLAGS_sanitize = -O1 -g3 -fsanitize=address,undefined -fno-omit-frame-pointer
LDFLAGS_sanitize  = -fsanitize=address,undefined

CXXFLAGS := $(COMMON) $(CXXFLAGS_$(MODE))
LDFLAGS  := $(LIBS) $(LDFLAGS_$(MODE))

# ---------- Folders ----------
SRC_DIR   := src
BUILD_DIR := build
PLOT_DIR  := scripts
OBJ_DIR   := $(BUILD_DIR)/obj
BIN_DIR   := $(BUILD_DIR)/bin

# ---------- Sources & objects ----------
SRCS_SRC   := $(wildcard $(SRC_DIR)/*.cpp)              # library/engine sources
OBJS_COMMON:= $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS_SRC))

# Program entry points (in project root)
ENTRY_FEW  := few_particles.cpp
ENTRY_N    := N_particles.cpp

OBJ_FEW    := $(OBJ_DIR)/$(notdir $(ENTRY_FEW:.cpp=.o))
OBJ_N      := $(OBJ_DIR)/$(notdir $(ENTRY_N:.cpp=.o))

DEPS       := $(OBJS_COMMON:.o=.d) $(OBJ_FEW:.o=.d) $(OBJ_N:.o=.d)

# ---------- Executables ----------
EXE_FEW := $(BIN_DIR)/few_particles.exe
EXE_N   := $(BIN_DIR)/N_particles.exe
PLOT_SCRIPT := $(PLOT_DIR)/plot_pos_vel.py

.PHONY: all clean run-few run-n debug release sanitize dirs plot

# Build both by default
all: $(EXE_FEW) $(EXE_N)

debug:    MODE=debug
debug:    clean all

release:  MODE=release
release:  clean all

sanitize: MODE=sanitize
sanitize: clean all

dirs:
	@mkdir -p $(OBJ_DIR) $(BIN_DIR)

# ---------- Compilation rules ----------
# Compile objects from src/
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | dirs
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile entry points in project root
$(OBJ_DIR)/%.o: %.cpp | dirs
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---------- Link rules ----------
$(EXE_FEW): $(OBJS_COMMON) $(OBJ_FEW) | dirs
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(EXE_N): $(OBJS_COMMON) $(OBJ_N) | dirs
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# ---------- Convenience targets ----------
run-few: $(EXE_FEW)
	$(EXE_FEW) $(ARGS)

run-n: $(EXE_N)
	$(EXE_N) $(ARGS)

plot:
	@echo "Generating plots"
	python3 $(PLOT_SCRIPT)

clean:
	@rm -rf $(BUILD_DIR)

-include $(DEPS)
