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
SRCS_SRC    := $(wildcard $(SRC_DIR)/*.cpp)                             # library/engine sources
OBJS_COMMON := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS_SRC))

# Program entry points (in project root)
ENTRY_FEW  := few_particles.cpp
ENTRY_N    := N_particles.cpp
ENTRY_S    := single_particle.cpp

OBJ_FEW    := $(OBJ_DIR)/$(notdir $(ENTRY_FEW:.cpp=.o))
OBJ_N      := $(OBJ_DIR)/$(notdir $(ENTRY_N:.cpp=.o))
OBJ_S      := $(OBJ_DIR)/$(notdir $(ENTRY_S:.cpp=.o))

DEPS       := $(OBJS_COMMON:.o=.d) $(OBJ_FEW:.o=.d) $(OBJ_N:.o=.d) $(OBJ_S:.o=.d)

# ---------- Executables ----------
EXE_FEW    := $(BIN_DIR)/few_particles.exe
EXE_N      := $(BIN_DIR)/N_particles.exe
EXE_S      := $(BIN_DIR)/single_particle.exe
PLOT_FEW   := $(PLOT_DIR)/plot_pos_vel.py
PLOT_SINGLE:= $(PLOT_DIR)/plot_single.py

# ---------- Phony targets ----------
.PHONY: all clean run-few run-n run-s debug release sanitize dirs plot help print-%

# ---------- Default ----------
all: $(EXE_FEW) $(EXE_N) ## Build the main executables (few & N) in the selected MODE

# ---------- Build modes ----------
debug:    MODE=debug
debug:    clean all ## Rebuild everything with debug flags

release:  MODE=release
release:  clean all ## Rebuild everything with release flags

sanitize: MODE=sanitize
sanitize: clean all ## Rebuild with Address/UBSan enabled

# ---------- Utility ----------
dirs: ## Create build directories if missing
	@mkdir -p $(OBJ_DIR) $(BIN_DIR)

# Print any variable: make print-VAR (e.g., make print-CXXFLAGS)
print-%: ## Print a Makefile variable (usage: make print-VARNAME)
	@printf "%s = %s\n" "$*" "$($*)"

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

$(EXE_S): $(OBJS_COMMON) $(OBJ_S) | dirs
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# ---------- Run helpers ----------
run-few: $(EXE_FEW) ## Run few_particles executable (pass ARGS="..." to forward CLI args)
	$(EXE_FEW) $(ARGS)

run-n: $(EXE_N) ## Run N_particles executable (pass ARGS="...")
	$(EXE_N) $(ARGS)

run-s: $(EXE_S) ## Run single_particle executable (pass ARGS="...")
	$(EXE_S) $(ARGS)

# ---------- Plotting ----------
plot: ## Generate plots using the Python scripts in scripts/
	@echo "Generating plots"
	python3 $(PLOT_FEW)
	python3 $(PLOT_SINGLE)

# ---------- Cleaning ----------
clean: ## Remove the entire build/ directory
	@rm -rf $(BUILD_DIR)

# ---------- Self-documenting help ----------
# Any target line with a double-hash comment (##) will be picked up by 'make help'.
help: ## Show this help message with available targets and current config
	@printf "Usage: make [TARGET] [MODE=debug|release|sanitize] [ARGS=\"...\"]\n\n"
	@printf "Current configuration:\n"
	@printf "  MODE     : %s\n" "$(MODE)"
	@printf "  CXX      : %s\n" "$(CXX)"
	@printf "  CXXFLAGS : %s\n" "$(CXXFLAGS)"
	@printf "  LDFLAGS  : %s\n\n" "$(LDFLAGS)"
	@awk 'BEGIN {FS":.*##"; printf "Available targets:\n"} \
	     /^[a-zA-Z0-9_\/\.\-%]+:.*##/ {printf "  %-18s %s\n", $$1, $$2}' $(MAKEFILE_LIST)
