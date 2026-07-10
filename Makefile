# ===========================================================================
#  McLiCS - Makefile
#  Targets:
#    make            -> GPU binary (mc_sim), requires nvcc
#    make CPU        -> CPU-only binary (mc_sim)
#    make clean      -> remove build artefacts and the binary
#    make clean-data -> remove simulation output (ic.csv, po.dat, director_field_*.csv)
# ===========================================================================
 PROGRAM := mc_sim
########### Gnu:
 COMPILER := g++ 
 GPUCOMP  := nvcc
 FLAGS    :=  -O3  -fopenmp -static
 GPUFLAGS    :=  -O3
 #-Wno-unused-result -Wno-format
 LIB := -lm  -lgsl -lgslcblas  -lgomp 
# LIB := -lm  -lgslcblas  -lgsl -lgomp 
USECUDA:=-D CUDA__="CUDA"
ifeq (,$(shell which nvcc))
	.DEFAULT_GOAL=CPU
	USECUDA:= 
endif
############ Intel 2 :
#COMPILER = icpc
#FLAGS= -ipo -O3  -no-prec-div -xAVX -simd -qopenmp -fp-model fast=2 -static
#LIB = -mkl -lgsl 


############ Intel Debug :

#COMPILER = icpc
#FLAGS= -O0 -g -static
#LIB = -mkl -lgsl 



 
CPPS := $(wildcard src/*.cpp)
CUDA := $(wildcard src/*.cu)
HEADER := $(wildcard include/*.h)
OBJS  := $(patsubst src/%.cpp,build/%.o,${CPPS})
CBJS  := $(patsubst src/%.cu,build/%.cuda.o,${CUDA})
DBGOBJS  := $(patsubst src/%.cpp,build/%dbg.o,${CPPS})
DBGCBJS  := $(patsubst src/%.cu,build/%dbg.cuda.o,${CUDA})

all:${PROGRAM}

${PROGRAM}: ${OBJS} ${CBJS}
	@${GPUCOMP} ${GPUFLAGS} $(filter-out  build/simulator.o,${OBJS}) ${CBJS} ${LIB} -o ${PROGRAM}
CPU:  $(filter-out  build/simulatorGPU.o,${OBJS})
	@${COMPILER} ${FLAGS} $(filter-out  build/simulatorGPU.o,${OBJS}) ${LIB} -o ${PROGRAM}
${OBJS}: build/%.o: src/%.cpp ${HEADER} | build
	${COMPILER} ${FLAGS} -c $< -o $@ ${USECUDA}
${CBJS}: build/%.cuda.o: src/%.cu  | build
	${GPUCOMP}  ${GPUFLAGS} -dc -c  $< -o $@ -D CUDA__="CUDA"
	
build:
	@mkdir build
debug:	${DBGOBJS} ${DBGCBJS}
	${GPUCOMP}  -O0 -g -Xcompiler -fopenmp  -lineinfo $(filter-out -O% -fast -static, ${GPUFLAGS}) $(filter-out  build/simulatordbg.o,${DBGOBJS}) ${DBGCBJS}  ${LIB} -o ${PROGRAM}_debug
	#@rm build/*dbg.o
${DBGOBJS}: build/%dbg.o: src/%.cpp  | build
	${COMPILER}  -O0 -g $(filter-out -O% -fast, ${FLAGS}) -c $< -o $@ ${USECUDA}
${DBGCBJS}: build/%dbg.cuda.o: src/%.cu | build
	${GPUCOMP}  -O0 -g $(filter-out -O% -fast, ${GPUFLAGS}) -lineinfo -dc -c $< -o $@ ${USECUDA}
$(patsubst %.o,%dbg.o,${OBJS}): ${HEADER}

renew: clean	${PROGRAM}

clean:
	@rm -f ${PROGRAM}
	@rm -fr build
	@rm -f *.o
	@echo "Clean done."

clean-data:
	@rm -f ic.csv po.dat
	@find . -maxdepth 1 -type f -name "director_field_*.csv" -delete
	@echo "Simulation output removed (ic.csv, po.dat, director_field_*.csv)."

.PHONY: all CPU debug build clean clean-data renew
