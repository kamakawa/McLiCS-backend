# Nome dos executáveis
PROGRAM_GPU := mc_sim_gpu
PROGRAM_CPU := mc_sim_cpu

# Diretórios
SRCDIR := rod/src
BUILDDIR := rod/build
INCDIR := rod/include

# Compilador e bibliotecas
COMPILER := g++
GPUCOMP := nvcc
FLAGS := -O3 -fopenmp -static
LIB := -lm -lgsl -lgslcblas -lgomp
HEADER_FLAGS := -I$(INCDIR)

# Arquivos de código-fonte
CPPS_ALL := $(wildcard $(SRCDIR)/*.cpp)
CUDA_ALL := $(wildcard $(SRCDIR)/*.cu)
HEADER_ALL := $(wildcard $(INCDIR)/*.h)

# Gerar a lista de arquivos de objeto (.o)
# Versão CPU
CPPS_CPU := $(filter-out $(SRCDIR)/simulatorGPU.cpp, ${CPPS_ALL})
OBJS_CPU := $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,${CPPS_CPU})

# Versão GPU
CPPS_GPU := $(filter-out $(SRCDIR)/simulator.cpp, ${CPPS_ALL})
OBJS_GPU := $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,${CPPS_GPU})
CBJS_GPU := $(patsubst $(SRCDIR)/%.cu,$(BUILDDIR)/%.cuda.o,${CUDA_ALL})

# Define o alvo padrão
.DEFAULT_GOAL := CPU

# Regra para compilar e vincular a versão de CPU
CPU: ${PROGRAM_CPU}
${PROGRAM_CPU}: ${OBJS_CPU}
	@echo "Vincular objetos para a versão de CPU..."
	@${COMPILER} ${FLAGS} ${OBJS_CPU} ${LIB} -o $@

# Regra para compilar e vincular a versão de GPU
GPU: ${PROGRAM_GPU}
${PROGRAM_GPU}: ${OBJS_GPU} ${CBJS_GPU}
	@echo "Vincular objetos para a versão de GPU..."
	@${GPUCOMP} ${FLAGS} ${OBJS_GPU} ${CBJS_GPU} ${LIB} -o $@

# Regra para compilar arquivos C++ em .o
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp ${HEADER_ALL} | $(BUILDDIR)
	@echo "Compilando arquivo C++ $<"
	@${COMPILER} ${FLAGS} ${HEADER_FLAGS} -c $< -o $@

# Regra para compilar arquivos CUDA em .o
$(BUILDDIR)/%.cuda.o: $(SRCDIR)/%.cu ${HEADER_ALL} | $(BUILDDIR)
	@echo "Compilando arquivo CUDA $<"
	@${GPUCOMP} ${HEADER_FLAGS} -dc -c $< -o $@

# Cria o diretório de build
$(BUILDDIR):
	@mkdir -p $(BUILDDIR)

# Limpa o projeto
clean:
	@echo "Limpando arquivos do projeto..."
	@rm -f ${PROGRAM_CPU} ${PROGRAM_GPU}
	@rm -fr ${BUILDDIR}

renew: clean CPU