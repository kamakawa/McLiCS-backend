# ===========================================================================
#  McLiCS — Makefile
#  Targets:
#    make          → GPU binary (mc_sim),  requires nvcc
#    make cpu      → CPU binary (mc_sim_cpu)
#    make all      → both binaries (if nvcc is available)
#    make debug    → GPU binary with debug symbols
#    make debug_cpu→ CPU binary with debug symbols
#    make clean    → remove build artefacts and binaries
# ===========================================================================

PROGRAM_GPU := mc_sim
PROGRAM_CPU := mc_sim_cpu

COMPILER := g++
GPUCOMP  := nvcc

# Detect architecture from the GPU at build time (falls back to sm_60)
GPU_ARCH := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null \
            | head -1 | tr -d '.' | sed 's/^/sm_/' || echo sm_60)

FLAGS    := -O3 -std=c++17 -fopenmp
GPUFLAGS := -O3 -std=c++17 -arch=$(GPU_ARCH)

LIB      := -lm -lgsl -lgslcblas -lgomp

# Detect whether nvcc is available
HAS_CUDA := $(shell which nvcc > /dev/null 2>&1 && echo yes)

# ---------------------------------------------------------------------------
# Source lists
# simulator.cpp is the single unified simulator (replaces simulatorGPU.cpp)
# ---------------------------------------------------------------------------
ALL_CPPS := $(wildcard src/*.cpp)

CUDA_SRC := $(wildcard src/*.cu)
CUDA_OBJ := $(patsubst src/%.cu,  build/%.cuda.o, $(CUDA_SRC))

CPU_OBJ  := $(patsubst src/%.cpp, build/cpu_%.o,  $(ALL_CPPS))
GPU_OBJ  := $(patsubst src/%.cpp, build/gpu_%.o,  $(ALL_CPPS))

HEADERS  := $(wildcard include/*.h) $(wildcard include/*.cuh)

# ---------------------------------------------------------------------------
# Default target
# ---------------------------------------------------------------------------
ifeq ($(HAS_CUDA),yes)
.DEFAULT_GOAL := gpu
else
.DEFAULT_GOAL := cpu
endif

# ---------------------------------------------------------------------------
# GPU build
# ---------------------------------------------------------------------------
ifeq ($(HAS_CUDA),yes)

gpu: build $(PROGRAM_GPU)

$(PROGRAM_GPU): $(GPU_OBJ) $(CUDA_OBJ)
	@echo "=== Linking GPU binary ($(GPU_ARCH)) ==="
	$(GPUCOMP) $(GPUFLAGS) $^ $(LIB) -o $@
	@echo "GPU binary ready: ./$(PROGRAM_GPU)"

build/gpu_%.o: src/%.cpp $(HEADERS) | build
	@echo "  [GPU-cpp] $<"
	$(COMPILER) $(FLAGS) -DCUDA__="CUDA" -c $< -o $@

build/%.cuda.o: src/%.cu $(HEADERS) | build
	@echo "  [CUDA]    $<"
	$(GPUCOMP) $(GPUFLAGS) -DCUDA__="CUDA" -dc -c $< -o $@

else

gpu:
	@echo "Error: nvcc not found. Cannot build GPU target."
	@exit 1

endif

# ---------------------------------------------------------------------------
# CPU build
# ---------------------------------------------------------------------------
cpu: build $(PROGRAM_CPU)

$(PROGRAM_CPU): $(CPU_OBJ)
	@echo "=== Linking CPU binary ==="
	$(COMPILER) $(FLAGS) $^ $(LIB) -o $@
	@echo "CPU binary ready: ./$(PROGRAM_CPU)"

build/cpu_%.o: src/%.cpp $(HEADERS) | build
	@echo "  [CPU-cpp] $<"
	$(COMPILER) $(FLAGS) -march=native -c $< -o $@

# ---------------------------------------------------------------------------
# Debug targets
# ---------------------------------------------------------------------------
debug: build
ifeq ($(HAS_CUDA),yes)
	@echo "=== Debug GPU build ==="
	$(GPUCOMP) -O0 -g -G -Xcompiler -fopenmp -lineinfo \
	    -DCUDA__="CUDA" $(GPU_OBJ) $(CUDA_OBJ) $(LIB) -o $(PROGRAM_GPU)_debug
else
	@echo "Error: nvcc not found."
	@exit 1
endif

debug_cpu: $(CPU_OBJ) | build
	@echo "=== Debug CPU build ==="
	$(COMPILER) -O0 -g $(FLAGS) $(ALL_CPPS) $(LIB) -o $(PROGRAM_CPU)_debug

# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------
build:
	@mkdir -p build

all: cpu
ifeq ($(HAS_CUDA),yes)
all: gpu
endif

clean:
	@rm -f $(PROGRAM_GPU) $(PROGRAM_CPU) $(PROGRAM_GPU)_debug $(PROGRAM_CPU)_debug
	@rm -rf build
	@echo "Clean done."

clean-data:
	@find . -type f -name "director_field_*.csv" -delete
	@echo "All director_field_*.csv files removed."
	
.PHONY: gpu cpu all debug debug_cpu clean clean-data
