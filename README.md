# 🧊 Simulador Monte Carlo 3D para Cristais Líquidos com Ancoramento Superficial

> Projeto de Iniciação Científica desenvolvido por **Eric Kamakawa**, sob orientação dos professores **Rafael Zola** e **Rodolfo Teixeira**, no curso de Engenharia da Computação da UTFPR – Apucarana.

---

## 🎯 Objetivo

Este projeto de Iniciação Científica visa dar continuidade e aprimorar o simulador computacional McLiCS (Monte Carlo Liquid Crystal Simulator), explorando o comportamento de sistemas de cristais líquidos nemáticos em três dimensões. O foco principal é a implementação e análise do **modelo de Lebwohl-Lasher** e a incorporação de diversas condições de ancoramento superficial para simulações em diferentes regimes físicos.

O projeto é uma aplicação direta de **Física Computacional**, combinando conceitos de física estatística e programação para resolver problemas complexos da matéria condensada.

---

## 🔬 Contexto Físico (Revisado e Expandido)

Cristais líquidos nemáticos representam um estado da matéria que exibe ordem orientacional sem ordem posicional de longo alcance, encontrando aplicações em tecnologias de ponta como telas LCD, sensores térmicos e materiais fotônicos.

O comportamento coletivo dessas moléculas é investigado utilizando o **Modelo de Lebwohl-Lasher (LL)**. Este é um modelo de rede (lattice) que simplifica cada molécula como um vetor tridimensional ($\mathbf{n}$) em um sítio de grade cúbica. A energia de interação entre moléculas vizinhas ($U_{ij}$) é descrita pela força de Van der Waals e é dada pelo potencial:

$$U_{ij}=-JP_2(\cos\theta_{ij})$$

onde $P_2(\cos\theta_{ij})$ é o segundo polinômio de Legendre e $\theta_{ij}$ é o ângulo entre as orientações dos vetores diretores vizinhos. Esta interação impõe um alinhamento paralelo (nemático) que é a força motriz para a transição de fase isotrópico-nemático.

O projeto também explora o papel crucial das **Superfícies (Ancoramento)**, que impõem restrições de contorno que afetam a orientação das moléculas. A energia de ancoramento é fundamental para a formação de texturas e defeitos topológicos, fenômenos que são críticos para a funcionalidade de dispositivos ópticos.

---

## 🧮 Metodologia Computacional (Revisada e Expandida)

A base do simulador é o **Método de Monte Carlo** no formalismo de **Metropolis**, uma técnica estocástica que simula o equilíbrio termodinâmico de sistemas. O processo envolve a perturbação da orientação de uma molécula e a aceitação ou rejeição dessa nova configuração com base na variação da energia total ($\Delta E$) e na temperatura ($T$), de acordo com a probabilidade de Boltzmann.

As simulações são executadas em um **lattice 3D** e focam na análise do **Parâmetro de Ordem Uniaxial ($\mathbf{S}$)**, que é a métrica fundamental para quantificar o grau de alinhamento do sistema.

A simulação é conduzida sob diversos regimes de evolução:

- **Evolução Térmica (Thermal Quench):** O sistema é resfriado gradualmente de uma temperatura inicial alta ($T_{i}$, estado isotrópico) até uma temperatura final baixa ($T_{f}$, estado nemático). O objetivo é mapear a **Transição de Fase** através do comportamento do parâmetro $S$ em função de $T$.
- **Análise em Temperatura Constante:** O sistema é mantido em uma temperatura específica por longos ciclos para garantir que o **equilíbrio térmico** seja atingido e propriedades estáticas sejam calculadas com precisão estatística.
- **Estudo de Campos Externos:** O potencial de interação é modificado para incluir termos de energia de **campo elétrico** (ou magnético), permitindo a análise da resposta eletro-óptica do material.
- **Condições de Ancoramento:** Implementação de diferentes modelos de energia de ancoramento (e.g., Rapine-Papoular, Homeotrópico, Forte), fundamentais para simular células de cristal líquido e estudar a formação de **domínios** e **defeitos**.

---

## 🏛️ Arquitetura do Software e POO

O projeto é estruturado em C++ moderno seguindo os princípios da **Programação Orientação a Objetos (POO)** para garantir modularidade e segurança:

- **Polimorfismo e Segurança:** Implementação de **destrutores virtuais** nas classes base (`Evolve`, `Geometry`) e uso consistente da palavra-chave `override`. Esta abordagem garante a correta desalocação de memória (evitando *memory leaks*) e a correta execução dos métodos nas classes derivadas, essenciais para o polimorfismo.
- **Coesão e Encapsulamento:** Uso de `structs` aninhadas para organizar logicamente os parâmetros da simulação (ex: `params->lattice.Nx`), e utilização de `namespaces` (ex: `namespace IO`) para isolar as funcionalidades e evitar a poluição do escopo global.
- **Herança:** O código é construído sobre uma hierarquia de classes (ex: `EvolveN` herda de `Evolve`), permitindo que novos métodos de evolução sejam adicionados facilmente.

---

## 🛠️ Tecnologias e Fluxo de Trabalho

### **Controle de Simulação (`input_parameters.txt`)**
O arquivo `input_parameters.txt` serve como a interface de controle do simulador. Ele permite que o pesquisador ajuste facilmente parâmetros críticos como a **temperatura inicial (Ti)**, a **temperatura final (Tf)** e o **passo de resfriamento (dT)** sem a necessidade de recompilar o código. Isso garante alta flexibilidade na exploração de diferentes cenários termodinâmicos.

### **Gerenciamento de Build (`Makefile`)**
O `Makefile` automatiza e gerencia todo o processo de **Computação Científica** do projeto. Ele é responsável por:
1.  Compilar os arquivos C++ e CUDA com `g++` e `nvcc`.
2.  Garantir a linkagem correta com bibliotecas externas essenciais (**GSL** e **OpenMP**).
3.  Permitir a construção separada das versões **CPU** (`make CPU`) e **GPU** (`make GPU`), otimizando o fluxo de trabalho.

### **Tecnologias Utilizadas**

- **Linguagem:** **C++** moderno (para performance).
- **Bibliotecas:** **GSL** (GNU Scientific Library - para números aleatórios e matemática) e **OpenMP/CUDA** (para paralelismo).
- **Análise de Dados:** **Python** com **Pandas** e **Matplotlib** (para análise e visualização dos resultados).

---

## 📚 Referências

- Tese de Doutorado – Eric Koudhi Omori (2018)
- Dissertação de Mestrado – Eric Koudhi Omori (2016)
- Landau & Binder – *A Guide to Monte Carlo Simulations in Statistical Physics*

---

## 📌 Futuras Etapas

- Implementação de CUDA para paralelismo em GPU
- Documentação
- Estudo de propriedades óticas e defeitos topológicos
- Escrita de artigo científico com base nos resultados

---

**Desenvolvido com ciência, código e propósito.**
