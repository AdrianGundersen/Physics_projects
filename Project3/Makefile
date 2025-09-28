# ---------- Compiler & flags ----------
CXX       ?= g++
CXXSTD     = -std=c++20
INC        = -Iinclude
LIBS       = -larmadillo

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
OBJ_DIR   := $(BUILD_DIR)/obj
BIN_DIR   := $(BUILD_DIR)/bin

# ---------- Sources & objects ----------
SRCS_SRC  := $(wildcard $(SRC_DIR)/*.cpp)
SRCS_MAIN := main.cpp
# flatten: src/foo.cpp → build/obj/foo.o
OBJS      := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS_SRC)) \
             $(OBJ_DIR)/$(SRCS_MAIN:.cpp=.o)
DEPS      := $(OBJS:.o=.d)

# ---------- Executable ----------
EXE := $(BIN_DIR)/project3.exe

.PHONY: all clean run debug release sanitize dirs

all: $(EXE)

debug:    MODE=debug
debug:    clean all

release:  MODE=release
release:  clean all

sanitize: MODE=sanitize
sanitize: clean all

dirs:
	@mkdir -p $(OBJ_DIR) $(BIN_DIR)

# Compile objects from src/
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | dirs
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile main.cpp from project root
$(OBJ_DIR)/main.o: main.cpp | dirs
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link final executable
$(EXE): $(OBJS) | dirs
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

run: $(EXE)
	$(EXE) $(ARGS)

clean:
	@rm -rf $(BUILD_DIR)

-include $(DEPS)
