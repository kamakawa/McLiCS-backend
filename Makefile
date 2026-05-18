PROGRAM_GPU := mc_sim
PROGRAM_CPU := mc_sim_cpu

########### Gnu:
COMPILER := g++
GPUCOMP  := nvcc
FLAGS    := -O3 -std=c++17 -fopenmp
GPUFLAGS := -O3 -std=c++17 -arch=sm_60  # Ajuste para sua GPU
LIB      := -lm -lgsl -lgslcblas -lgomp

# Detecta se nvcc está disponível
HAS_CUDA := $(shell which nvcc > /dev/null 2>&1 && echo yes)

# Exclui simulatorGPU.cpp e simulator.cpp (usamos apenas um)
ALL_CPPS := $(wildcard src/*.cpp)

# Para CPU: use apenas simulator.cpp (versão híbrida), exclua simulatorGPU.cpp
CPPS_FOR_CPU := $(filter-out src/simulatorGPU.cpp,${ALL_CPPS})

# Para GPU: use apenas simulatorGPU.cpp, exclua simulator.cpp
CPPS_FOR_GPU := $(filter-out src/simulator.cpp,${ALL_CPPS})

HEADER   := $(wildcard include/*.h)

# Objetos para CPU
OBJS_CPU := $(patsubst src/%.cpp,build/cpu_%.o,${CPPS_FOR_CPU})

# Objetos para GPU
OBJS_GPU := $(patsubst src/%.cpp,build/gpu_%.o,${CPPS_FOR_GPU})

# Arquivos CUDA
CUDA     := $(wildcard src/*.cu)
CBJS     := $(patsubst src/%.cu,build/%.cuda.o,${CUDA})

all: ${PROGRAM_GPU}

# -------------------- GPU build (usa simulatorGPU.cpp) --------------------
ifeq ($(HAS_CUDA),yes)
    ${PROGRAM_GPU}: ${OBJS_GPU} ${CBJS} | build
	@echo "========== Building GPU version =========="
	@${GPUCOMP} ${GPUFLAGS} ${OBJS_GPU} ${CBJS} ${LIB} -o ${PROGRAM_GPU}
	@echo "GPU executable: ./${PROGRAM_GPU}"
	
    build/gpu_%.o: src/%.cpp ${HEADER} | build
	@echo "Compiling $< for GPU build..."
	${COMPILER} ${FLAGS} -c $< -o $@ -D CUDA__="CUDA"
	
    build/%.cuda.o: src/%.cu | build
	@echo "Compiling CUDA $<..."
	${GPUCOMP} ${GPUFLAGS} -dc -c $< -o $@ -D CUDA__="CUDA"
else
    ${PROGRAM_GPU}:
	@echo "CUDA not available. Cannot build GPU version."
	@exit 1
endif

# -------------------- CPU build (usa simulator.cpp com CPU paths) --------------------
CPU: ${OBJS_CPU} | build
	@echo "========== Building CPU version =========="
	@${COMPILER} ${FLAGS} ${OBJS_CPU} ${LIB} -o ${PROGRAM_CPU}
	@echo "CPU executable: ./${PROGRAM_CPU}"

build/cpu_%.o: src/%.cpp ${HEADER} | build
	@echo "Compiling $< for CPU build..."
	${COMPILER} ${FLAGS} -c $< -o $@

build:
	@mkdir -p build

# -------------------- Debug --------------------
debug: ${OBJS_GPU} ${CBJS} | build
	@echo "========== Building DEBUG version =========="
	${GPUCOMP} -O0 -g -G -Xcompiler -fopenmp -lineinfo \
	${GPUFLAGS} ${OBJS_GPU} ${CBJS} ${LIB} -o ${PROGRAM_GPU}_debug

# -------------------- Clean --------------------
clean:
	@rm -f ${PROGRAM_GPU} ${PROGRAM_CPU} ${PROGRAM_GPU}_debug
	@rm -fr build
	@rm -f *.o

renew: clean CPU

.PHONY: all CPU debug clean renew