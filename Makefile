PROGRAM_GPU := mc_sim
PROGRAM_CPU := mc_sim_cpu

OUT_DIR := dist/linux-x64

########### Gnu:
COMPILER := g++
GPUCOMP  := nvcc
FLAGS    := -O3 -std=c++17 -fopenmp -static
GPUFLAGS := -O3 -std=c++17
LIB      := -lm -lgsl -lgslcblas -lgomp

USECUDA := -D CUDA__="CUDA"
ifeq (,$(shell which nvcc))
    .DEFAULT_GOAL=CPU
    USECUDA :=
endif

# simulatorGPU.cpp was merged into simulator.cpp — exclude it from all builds
ALL_CPPS := $(wildcard src/*.cpp)
CPPS     := $(filter-out src/simulatorGPU.cpp,${ALL_CPPS})
CUDA     := $(wildcard src/*.cu)
HEADER   := $(wildcard include/*.h)

OBJS    := $(patsubst src/%.cpp,build/%.o,${CPPS})
CBJS    := $(patsubst src/%.cu,build/%.cuda.o,${CUDA})
DBGOBJS := $(patsubst src/%.cpp,build/%dbg.o,${CPPS})
DBGCBJS := $(patsubst src/%.cu,build/%dbg.cuda.o,${CUDA})

all: ${PROGRAM_GPU}

# ---------------- GPU build ----------------
${PROGRAM_GPU}: ${OBJS} ${CBJS} | ${OUT_DIR}
	@${GPUCOMP} ${GPUFLAGS} ${OBJS} ${CBJS} ${LIB} -o ${OUT_DIR}/${PROGRAM_GPU}

# ---------------- CPU build ----------------
CPU: ${OBJS} | ${OUT_DIR}
	@${COMPILER} ${FLAGS} ${OBJS} ${LIB} -o ${OUT_DIR}/${PROGRAM_CPU}

# -------------- object rules --------------
${OBJS}: build/%.o: src/%.cpp ${HEADER} | build
	${COMPILER} ${FLAGS} -c $< -o $@ ${USECUDA}

${CBJS}: build/%.cuda.o: src/%.cu | build
	${GPUCOMP} ${GPUFLAGS} -dc -c $< -o $@ -D CUDA__="CUDA"

build:
	@mkdir -p build

${OUT_DIR}:
	@mkdir -p ${OUT_DIR}

# ---------------- debug ----------------
debug: ${DBGOBJS} ${DBGCBJS} | ${OUT_DIR}
	${GPUCOMP} -O0 -g -Xcompiler -fopenmp -lineinfo \
	$(filter-out -O% -fast -static, ${GPUFLAGS}) \
	${DBGOBJS} ${DBGCBJS} ${LIB} \
	-o ${OUT_DIR}/${PROGRAM_GPU}_debug

${DBGOBJS}: build/%dbg.o: src/%.cpp | build
	${COMPILER} -O0 -g $(filter-out -O% -fast, ${FLAGS}) -c $< -o $@ ${USECUDA}

${DBGCBJS}: build/%dbg.cuda.o: src/%.cu | build
	${GPUCOMP} -O0 -g $(filter-out -O% -fast, ${GPUFLAGS}) -lineinfo -dc -c $< -o $@ ${USECUDA}

$(patsubst %.o,%dbg.o,${OBJS}): ${HEADER}

renew: clean all

clean:
	@rm -f ${OUT_DIR}/${PROGRAM_GPU} ${OUT_DIR}/${PROGRAM_CPU} ${OUT_DIR}/${PROGRAM_GPU}_debug
	@rm -fr build
	@rm -f *.o

clean-data:
	@echo "Removing output data files: director_field_*.csv"
	@rm -f director_field_*.csv
