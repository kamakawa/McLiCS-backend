# Nome dos executáveis
PROGRAM_CPU := mc_sim_cpu

# Diretórios
SRCDIR := rod/src
BUILDDIR := rod/build
INCDIR := rod/include

# Compilador e bibliotecas
COMPILER := g++
FLAGS := -O3 -fopenmp -static
LIB := -lm -lgsl -lgslcblas -lgomp
HEADER_FLAGS := -I$(INCDIR)

# Arquivos de código-fonte
# Lista TODOS os arquivos .cpp no diretório src, excluindo o simulatorGPU (se existir)
CPPS_ALL := $(shell find $(SRCDIR) -name '*.cpp')
CPPS_CPU := $(filter-out $(SRCDIR)/simulatorGPU.cpp, ${CPPS_ALL})

# Gerar a lista de arquivos de objeto (.o)
# Transforma cada caminho de .cpp em um caminho de .o no diretório BUILDDIR
OBJS_CPU := $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,${CPPS_CPU})

# Lista de todos os headers para dependência de compilação
HEADER_ALL := $(wildcard $(INCDIR)/*.h)

# Define o alvo padrão
.DEFAULT_GOAL := CPU

# ----------------------------------------------------------------------
# Regras de Build
# ----------------------------------------------------------------------

# Regra para compilar e vincular a versão de CPU (PROGRAM_CPU)
CPU: ${PROGRAM_CPU}

${PROGRAM_CPU}: ${OBJS_CPU}
	@echo "Vincular objetos para a versão de CPU..."
	@${COMPILER} ${FLAGS} ${OBJS_CPU} ${LIB} -o $@

# Regra para compilar arquivos C++ em .o
# O pipe '| $(BUILDDIR)' garante que o diretório exista antes de compilar
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp ${HEADER_ALL} | $(BUILDDIR)
	@echo "Compilando arquivo C++ $<"
	@${COMPILER} ${FLAGS} ${HEADER_FLAGS} -c $< -o $@

# Cria o diretório de build
$(BUILDDIR):
	@mkdir -p $(BUILDDIR)

# ----------------------------------------------------------------------
# Regras de Limpeza
# ----------------------------------------------------------------------

# Limpa o projeto
clean:
	@echo "Limpando arquivos do projeto..."
	@rm -f ${PROGRAM_CPU}
	@rm -fr ${BUILDDIR}

# Alvo para limpar e recompilar
renew: clean CPU

# Fim do Makefile