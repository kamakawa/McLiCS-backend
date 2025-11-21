PROGRAM := mc_sim
########### Gnu:
COMPILER := g++ 
GPUCOMP  := nvcc
FLAGS    :=  -O3  -fopenmp -static
GPUFLAGS :=  -O3
LIB := -lm  -lgsl -lgslcblas  -lgomp 
USECUDA:=-D CUDA__="CUDA"
ifeq (,$(shell which nvcc))
	.DEFAULT_GOAL=CPU
	USECUDA:= 
endif

# ========== NOVOS ARQUIVOS ==========
CPPS := $(wildcard src/*.cpp)
CUDA := $(wildcard src/*.cu)
HEADER := $(wildcard include/*.h)

# REMOVER arquivos antigos de evolve
OLD_EVOLVE := src/evolve_thermal.cpp src/evolve_step.cpp src/evolve_quench.cpp src/evolve_electric.cpp src/evolve.cpp
CPPS := $(filter-out $(OLD_EVOLVE),$(CPPS))

OBJS  := $(patsubst src/%.cpp,build/%.o,${CPPS})
CBJS  := $(patsubst src/%.cu,build/%.cuda.o,${CUDA})
DBGOBJS  := $(patsubst src/%.cpp,build/%dbg.o,${CPPS})
DBGCBJS  := $(patsubst src/%.cu,build/%dbg.cuda.o,${CUDA})

all:${PROGRAM}

${PROGRAM}: ${OBJS} ${CBJS}
	@${GPUCOMP} ${GPUFLAGS} $(filter-out build/simulator.o,${OBJS}) ${CBJS} ${LIB} -o ${PROGRAM}

CPU: $(filter-out build/simulatorGPU.o,${OBJS})
	@${COMPILER} ${FLAGS} $(filter-out build/simulatorGPU.o,${OBJS}) ${LIB} -o ${PROGRAM}

${OBJS}: build/%.o: src/%.cpp ${HEADER} | build
	${COMPILER} ${FLAGS} -c $< -o $@ ${USECUDA}

${CBJS}: build/%.cuda.o: src/%.cu | build
	${GPUCOMP} ${GPUFLAGS} -dc -c $< -o $@ -D CUDA__="CUDA"
	
build:
	@mkdir build

debug: ${DBGOBJS} ${DBGCBJS}
	${GPUCOMP} -O0 -g -Xcompiler -fopenmp -lineinfo $(filter-out -O% -fast -static, ${GPUFLAGS}) $(filter-out build/simulatordbg.o,${DBGOBJS}) ${DBGCBJS} ${LIB} -o ${PROGRAM}_debug

${DBGOBJS}: build/%dbg.o: src/%.cpp | build
	${COMPILER} -O0 -g $(filter-out -O% -fast, ${FLAGS}) -c $< -o $@ ${USECUDA}

${DBGCBJS}: build/%dbg.cuda.o: src/%.cu | build
	${GPUCOMP} -O0 -g $(filter-out -O% -fast, ${GPUFLAGS}) -lineinfo -dc -c $< -o $@ ${USECUDA}

$(patsubst %.o,%dbg.o,${OBJS}): ${HEADER}

renew: clean ${PROGRAM}

clean:
	@rm -f ${PROGRAM}
	@rm -fr build
	@rm -f *.o

# ========== NOVA REGRA: Limpeza específica dos arquivos antigos ==========
clean-old-evolve:
	@echo "Removendo arquivos antigos de evolve..."
	@rm -f src/evolve_thermal.cpp src/evolve_step.cpp src/evolve_quench.cpp src/evolve_electric.cpp src/evolve.cpp
	@rm -f include/evolve.h
	@echo "Arquivos antigos removidos!"

# ========== VERIFICAÇÃO DE ARQUIVOS ==========
check-files:
	@echo "=== Arquivos CPP encontrados ==="
	@for file in ${CPPS}; do echo "$$file"; done
	@echo "=== Headers encontrados ==="  
	@for file in ${HEADER}; do echo "$$file"; done
	@echo "=== Objetos a serem gerados ==="
	@for file in ${OBJS}; do echo "$$file"; done