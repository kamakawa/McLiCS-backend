PROGRAM   := mc_sim
CXX       := g++
NVCC      := /usr/bin/nvcc

# Flags base (CPU puro)
CXXFLAGS      := -O3 -fopenmp -Iinclude
CXXFLAGS_CUDA := $(CXXFLAGS) -DUSE_CUDA

# Flags CUDA
NVCCFLAGS := -O3 -lineinfo -Iinclude -DUSE_CUDA

LIBS      := -lm -lgsl -lgslcblas -lgomp

HAVE_CUDA := $(shell command -v $(NVCC) >/dev/null 2>&1 && echo 1 || echo 0)

# fontes
CPPS  := $(wildcard src/*.cpp)
CUSRC := $(wildcard src/cuda/*.cu)

# headers (RECURSIVO, inclui include/cuda/*.cuh)
HEADERS := $(shell find include -type f \( -name "*.h" -o -name "*.hpp" -o -name "*.cuh" \))

# objetos
OBJS    := $(patsubst src/%.cpp,build/%.o,$(CPPS))
CUOBJS  := $(patsubst src/cuda/%.cu,build/cuda/%.cu.o,$(CUSRC))

ifeq ($(HAVE_CUDA),1)
.DEFAULT_GOAL := all
else
.DEFAULT_GOAL := cpu
endif

all: $(PROGRAM)
cpu: $(PROGRAM)_cpu

# ---------------- link ----------------
# GPU build (compila .cpp com USE_CUDA + linka .cu)
$(PROGRAM): CXXFLAGS := $(CXXFLAGS_CUDA)
$(PROGRAM): $(OBJS) $(CUOBJS) | build build/cuda
	$(NVCC) -O3 -lineinfo $(CUOBJS) $(OBJS) $(LIBS) -o $@

# CPU build (NÃO define USE_CUDA e não linka objetos CUDA)
$(PROGRAM)_cpu: CXXFLAGS := $(CXXFLAGS)
$(PROGRAM)_cpu: $(OBJS) | build
	$(CXX) $(CXXFLAGS) $(OBJS) $(LIBS) -o $@

# ---------------- compile cpp ----------------
build/%.o: src/%.cpp $(HEADERS) | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---------------- compile cu ----------------
build/cuda/%.cu.o: src/cuda/%.cu $(HEADERS) | build build/cuda
	$(NVCC) $(NVCCFLAGS) -dc -c $< -o $@

build:
	@mkdir -p build
build/cuda:
	@mkdir -p build/cuda

clean:
	@rm -rf build $(PROGRAM) $(PROGRAM)_cpu