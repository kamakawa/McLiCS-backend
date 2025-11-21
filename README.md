# 🧊 Simulador Monte Carlo 3D para Cristais Líquidos com Ancoramento Superficial  
### **Fase 2 da Iniciação Científica – Continuidade do Projeto McLiCS**

**Autor:** Eric Kamakawa  
**Orientadores:** Prof. Rafael Zola e Prof. Rodolfo Teixeira  
**Curso:** Engenharia da Computação – UTFPR, Apucarana  

---

## 🎯 Objetivo da Fase 2

Esta segunda fase da Iniciação Científica teve como objetivo **dar continuidade ao trabalho de Eric Koudhi Omori**, expandindo, corrigindo e reorganizando o simulador computacional **McLiCS (Monte Carlo Liquid Crystal Simulator)**, implementado em **C++ e CUDA**.

O foco central é simular o comportamento orientacional de **cristais líquidos nemáticos 3D**, utilizando o **Modelo de Lebwohl–Lasher** e diferentes condições de **ancoramento superficial**.  
Além disso, foram feitas melhorias estruturais no código, especialmente na arquitetura de classes, segurança de memória e modularização.

---

## 🔬 Contexto Físico

Cristais líquidos nemáticos formam um estado da matéria intermediário entre sólido e líquido, caracterizado pela **ordem orientacional** das moléculas.  
Para estudar esse sistema, o simulador utiliza o **Modelo de Lebwohl–Lasher (LL)**, onde cada molécula é representada por um vetor unitário $\mathbf{n}$ em uma rede cúbica 3D.

A energia de interação entre vizinhos é dada por:

\[
U_{ij} = -J P_2(\cos\theta_{ij}),
\]

onde:
- $P_2$ é o segundo polinômio de Legendre,
- $\theta_{ij}$ é o ângulo entre diretores vizinhos.

O simulador também implementa **superfícies com ancoramento**, fundamentais para reproduzir:
- alinhamento homeotrópico,
- alinhamento planar,
- variações de energia superficial,
- formação de defeitos topológicos,
- texturas ópticas relevantes para dispositivos LCD.

---

## 🧮 Metodologia Computacional

A simulação usa o **Método de Monte Carlo com o algoritmo de Metropolis**, em que cada passo envolve:

1. Selecionar um sítio aleatório da rede.  
2. Propor uma rotação aleatória no vetor diretor.  
3. Calcular a variação de energia $\Delta E$.  
4. Aceitar ou rejeitar com base no fator de Boltzmann:

\[
P = e^{-\Delta E/T}.
\]

O simulador trabalha em múltiplos regimes:

### 🔥 **1. Evolução Térmica (Thermal Quench)**
O sistema é resfriado de $T_i$ (isotrópico) até $T_f$ (nemático).  
Permite identificar a **transição de fase isotrópico–nemático**.

### 🧊 **2. Simulação em Temperatura Constante**
Modo de equilíbrio:  
o sistema é mantido em uma temperatura fixa para análise de propriedades estáticas.

### ⚡ **3. Campos Externos**
Termos adicionais no Hamiltoniano permitem simular:
- efeitos eletro-ópticos,
- realinhamento induzido por campo.

### 🧱 **4. Condições de Ancoramento Superficial**
A superfície impõe restrições à orientação das moléculas.  
O projeto implementa diferentes modelos de energia de ancoramento.

---

## 🏛️ Arquitetura e Organização em C++

Durante esta fase, o código recebeu melhorias importantes:

### ✔ **POO bem estruturada**
- Uso correto de *herança* entre classes como `Evolve`, `EvolveN` e `Geometry`.  
- Implementação de **destrutores virtuais**, fundamental para evitar *memory leaks*.  
- Métodos `override` para garantir segurança no polimorfismo.

### ✔ **Alta coesão e modularidade**
- Parâmetros organizados em `structs` dentro de namespaces.  
- Separação entre lógica física, geometria da rede e rotinas de evolução.  
- Criação de classes mais claras e expansíveis.

### ✔ **Namespaces e organização**
- Namespace `IO` para leitura/escrita.  
- Namespace `Utils` para funções auxiliares.

Essa abordagem tornou o simulador mais limpo, seguro e robusto para futuras expansões — inclusive para a implementação completa da GPU.

---

## 🛠️ Tecnologias Utilizadas

### **Linguagens**
- C++ (núcleo da simulação)
- CUDA (versão paralela em GPU – em avanço)
- Python (pós-processamento)

### **Bibliotecas**
- **GSL** – números aleatórios, matemática especializada  
- **OpenMP** – paralelismo CPU  
- **CUDA** – paralelismo massivo em GPU  
- **Pandas / Matplotlib** – análise e gráficos dos resultados  

---

## 🧩 Fluxo de Execução do Projeto

### 🔧 **1. Configuração pelo arquivo `input_parameters.txt`**
Esse arquivo controla toda a execução da simulação, incluindo:

- Dimensões da rede  
- Número de passos  
- Temperaturas (Ti, Tf, dT)  
- Intensidade de campo externo  
- Tipo de ancoramento  
- Arquivos de saída  

Permite fazer múltiplos experimentos sem recompilar o código.

---

### ⚙️ **2. Compilação via Makefile**
O Makefile gerencia todo o processo de build:

Com CPU:
```
make CPU
```

Com GPU:
```
make GPU
```

Executar:
```
./mclics
```

---

## 📈 Saídas Geradas

O simulador produz arquivos contendo:

- Parâmetro de ordem $S$  
- Energia média  
- Configurações orientacionais  
- Histogramas  
- Mapas 3D  

Esses dados são analisados posteriormente em Python.

---

## 📚 Referências Importantes

- Eric Koudhi Omori – Dissertação (2016) e Tese (2018)  
- Landau & Binder – *A Guide to Monte Carlo Simulations in Statistical Physics*  
- M. P. Allen & D. J. Tildesley – *Computer Simulation of Liquids*  

---

## 🚀 Próximos Passos

- Implementação completa da versão CUDA  
- Identificação e simulação de defeitos topológicos  
- Estudo de propriedades eletro-ópticas  
- Preparação de artigo científico  

---

**Desenvolvido com propósito, ciência e computação de alto desempenho.**  

