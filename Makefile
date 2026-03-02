# ============================================================
#                MClist - Monte Carlo Simulator
# ============================================================

PROGRAM := mc_sim

CXX     := g++
NVCC    := /usr/bin/nvcc

CXXFLAGS  := -O3 -fopenmp -Iinclude
NVCCFLAGS := -O3 -lineinfo -Iinclude -DUSE_CUDA
LIBS      := -lm -lgsl -lgslcblas -lgomp

# ------------------------------------------------------------
# Source files
# ------------------------------------------------------------
CPPS  := $(wildcard src/*.cpp)
CUSRC := $(wildcard src/cuda/*.cu)

HEADERS := $(shell find include -type f \( -name "*.h" -o -name "*.hpp" -o -name "*.cuh" \))

# ------------------------------------------------------------
# Objects
# ------------------------------------------------------------

# CPU objects (sem USE_CUDA)
OBJS_CPU := $(patsubst src/%.cpp,build_cpu/%.o,$(CPPS))

# GPU objects (com USE_CUDA)
OBJS_GPU := $(patsubst src/%.cpp,build_gpu/%.o,$(CPPS))
CUOBJS   := $(patsubst src/cuda/%.cu,build_gpu/cuda/%.cu.o,$(CUSRC))

# ------------------------------------------------------------
# Generated simulation files
# ------------------------------------------------------------
GENERATED_FILES := director_field_T_*.csv po.dat ic.csv

# ------------------------------------------------------------
# Default target
# ------------------------------------------------------------
.DEFAULT_GOAL := all

all: $(PROGRAM)
cpu: $(PROGRAM)_cpu

# ------------------------------------------------------------
# Link
# ------------------------------------------------------------

$(PROGRAM): $(OBJS_GPU) $(CUOBJS) | build_gpu build_gpu/cuda
	@echo "🔧 Linking GPU executable..."
	$(NVCC) -O3 -lineinfo $(CUOBJS) $(OBJS_GPU) $(LIBS) -o $@
	@echo "✅ GPU build complete: ./$(PROGRAM)"

$(PROGRAM)_cpu: $(OBJS_CPU) | build_cpu
	@echo "🔧 Linking CPU executable..."
	$(CXX) $(CXXFLAGS) $(OBJS_CPU) $(LIBS) -o $@
	@echo "✅ CPU build complete: ./$(PROGRAM)_cpu"

# ------------------------------------------------------------
# Compile C++ (CPU)
# ------------------------------------------------------------

build_cpu/%.o: src/%.cpp $(HEADERS) | build_cpu
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ------------------------------------------------------------
# Compile C++ (GPU mode)
# ------------------------------------------------------------

build_gpu/%.o: src/%.cpp $(HEADERS) | build_gpu
	$(CXX) $(CXXFLAGS) -DUSE_CUDA -c $< -o $@

# ------------------------------------------------------------
# Compile CUDA
# ------------------------------------------------------------

build_gpu/cuda/%.cu.o: src/cuda/%.cu $(HEADERS) | build_gpu build_gpu/cuda
	$(NVCC) $(NVCCFLAGS) -dc -c $< -o $@

# ------------------------------------------------------------
# Run targets
# ------------------------------------------------------------

run: $(PROGRAM)
	@echo "🚀 Running GPU simulation..."
	./$(PROGRAM) input_parameters.txt

run-cpu: $(PROGRAM)_cpu
	@echo "🚀 Running CPU simulation..."
	./$(PROGRAM)_cpu input_parameters.txt

# ------------------------------------------------------------
# Clean targets
# ------------------------------------------------------------

clean:
	@echo ""
	@echo "🔧 Cleaning build files..."
	@rm -rf build_cpu build_gpu $(PROGRAM) $(PROGRAM)_cpu
	@echo "✅ Build cleaned."
	@echo ""

clean-data:
	@echo ""
	@echo "🗑 Removing simulation output files..."
	@rm -f $(GENERATED_FILES)
	@echo "✅ Simulation data removed."
	@echo ""

purge: clean clean-data
	@echo ""
	@echo "🔥 Full cleanup completed."
	@echo ""

# ------------------------------------------------------------
# Create directories
# ------------------------------------------------------------

build_cpu:
	@mkdir -p build_cpu

build_gpu:
	@mkdir -p build_gpu

build_gpu/cuda:
	@mkdir -p build_gpu/cuda