<div align="center">

# McLiCS
### Monte Carlo Liquid Crystal Simulator

*Simulador Monte Carlo de cristais líquidos nemáticos, com núcleo em C++/CUDA.*

![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C?logo=cplusplus&logoColor=white)
![CUDA](https://img.shields.io/badge/CUDA-11%2B-76B900?logo=nvidia&logoColor=white)
![OpenMP](https://img.shields.io/badge/OpenMP-parallel-blue)
![GSL](https://img.shields.io/badge/GSL-2.x-orange)
![status](https://img.shields.io/badge/status-Iniciação%20Científica-informational)

**Autor:** Eric Kamakawa &nbsp;·&nbsp; **Orientadores:** Prof. Dr. Rafael Zola, Prof. Dr. Rodolfo Teixeira

**Curso:** Engenharia da Computação — UTFPR

</div>

---

## Sumário

- [Visão Geral](#visão-geral)
- [Dependências](#dependências)
- [Compilação](#compilação)
- [Uso](#uso)
- [Arquivo de Parâmetros](#arquivo-de-parâmetros)
- [Referência Completa de Parâmetros](#referência-completa-de-parâmetros)
- [Geometrias](#geometrias)
- [Modos de Evolução](#modos-de-evolução)
- [Potenciais](#potenciais)
- [Condições de Contorno](#condições-de-contorno)
- [Ancoragem de Superfície](#ancoragem-de-superfície)
- [Campo Elétrico](#campo-elétrico)
- [Saídas](#saídas)
- [Estrutura do Projeto](#estrutura-do-projeto)
- [Observações](#observações)

---

## Visão Geral

O **McLiCS** simula o comportamento de moléculas de cristal líquido nemático em uma rede cúbica tridimensional, utilizando o método de Monte Carlo com critério de Metropolis. A cada passo, o simulador calcula a energia de interação entre vizinhos próximos de cada sítio da rede e aceita ou rejeita perturbações angulares das moléculas.

Ao final de cada temperatura (ou cada passo/intensidade de campo, dependendo do modo), são calculados e salvos o **parâmetro de ordem escalar `S`** e a **energia média `E`** do sistema — além de um snapshot completo do campo diretor.

O núcleo de cálculo é escrito em **C++/CUDA** e compila em dois sabores:

| Build | Comando | Requer NVCC? | Suporta modos `*GPU`? |
|---|---|:---:|:---:|
| Completo (CPU + GPU) | `make` | Sim | Sim |
| Somente CPU | `make CPU` | Não | Não |

Em ambos os casos, o binário gerado é `./mc_sim`, e o paralelismo em CPU é feito via OpenMP.

---

## Dependências

| Dependência | Versão mínima | Necessário para |
|---|---|---|
| GCC | 9+ | Compilação (CPU e GPU) |
| NVCC (CUDA Toolkit) | 11+ | Compilação com suporte a GPU (`make`) |
| OpenMP | — | Paralelismo em CPU |
| GSL (GNU Scientific Library) | 2.x | Geração de números aleatórios e álgebra linear (autovalores do tensor Q) |

Instalação das dependências no Ubuntu/Debian:

```bash
sudo apt install g++ libgsl-dev libomp-dev
# Para suporte a GPU, instale o CUDA Toolkit compatível com seu driver:
# https://developer.nvidia.com/cuda-downloads
```

---

## Compilação

```bash
make            # binário completo (CPU + GPU), requer nvcc
make CPU        # binário somente CPU, não requer CUDA Toolkit
make clean      # remove binário e artefatos de build
make clean-data # remove saídas de simulações anteriores (ic.csv, po.dat, director_field_*.csv)
```

> O binário `mc_sim` gerado por `make` **também roda simulações em CPU normalmente** — o dispositivo usado em cada execução é decidido pelo parâmetro `evol` do arquivo de entrada (veja [Modos de Evolução](#modos-de-evolução)), não por uma flag de linha de comando.

---

## Uso

```bash
./mc_sim <arquivo_de_parâmetros>
```

**Exemplo:**

```bash
./mc_sim params.txt
```

Se nenhum arquivo for informado, o simulador roda com os valores padrão internos (ver tabelas abaixo).

---

## Arquivo de Parâmetros

O arquivo de parâmetros é um arquivo de texto com um parâmetro por linha, no formato `chave valor`. Linhas começando com `#` são comentários e são ignoradas. Os nomes dos parâmetros **não são case-sensitive**.

Exemplo completo:

```
# Rede
Nx          30
Ny          30
Nz          30
Nk          1                # obrigatório ser 1 (veja Observações)

# Geometria: bulk | slab | sphere | custom
geometry    bulk

# Condições de contorno: free | periodic
xbound      periodic
ybound      periodic
zbound      periodic

# Potencial: ll | ghrl | pear
potential   ll
A           1.0

# Modo de evolução: thermal | step | quench | electric (ou com sufixo GPU)
evol        thermal

# Temperatura inicial, final e passo
Ti          1.2
Tf          0.0
dT         -0.02

# Passos Monte Carlo de termalização e medição
MCT         500
MCS         200

# Condição inicial: random | homogeneous | ic_file | cholesteric
ic          random
```

---

## Referência Completa de Parâmetros

#### Rede

| Parâmetro | Valores | Padrão | Descrição |
|---|---|---|---|
| `Nx`, `Ny`, `Nz` | inteiro | 16 | Tamanho da rede em cada eixo |
| `Nk` (`neighbour_kind`) | `1` | 1 | Ordem de vizinhança considerada. **Só `1` é suportado** — outro valor encerra o programa com erro |
| `geometry` | `bulk` / `slab` / `sphere` / `custom` | `bulk` | Geometria da simulação |
| `xbound`, `ybound`, `zbound` | `free` / `periodic` | `free` | Condição de contorno por eixo |

#### Potencial e Constantes Elásticas

| Parâmetro | Valores | Padrão | Descrição |
|---|---|---|---|
| `potential` | `ll` / `ghrl` / `pear` | `ll` | Potencial de interação |
| `A` | float | 1.0 | Parâmetro do potencial LL |
| `B1`, `B2` | float | 0.04, 0.4 | Termos adicionais do potencial LL |
| `C` | float | 0.3 | Termo de acoplamento splay do potencial LL |
| `k11`, `k22`, `k33` | float | 1.0 | Constantes elásticas de Frank — splay, twist e bend (GHRL) |
| `p0` | float | 0 | Passo colestérico (GHRL; `0` = nemático puro) |

#### Evolução e Monte Carlo

| Parâmetro | Aliases | Padrão | Descrição |
|---|---|---|---|
| `evol` | `mode` | `thermal` | Modo de evolução — `thermal`, `step`, `quench`, `electric`, ou com sufixo `GPU` (ex. `thermalGPU`) para rodar em placa de vídeo |
| `Ti`, `Tf`, `dT` | — | 1.2 / 0.0 / -0.02 | Temperatura inicial, final e passo |
| `MCT` | — | 100 | Passos de termalização por temperatura |
| `MCS` | — | 100 | Passos de medição por temperatura |
| `fn` | — | 1 | Número de blocos independentes (modos `step`/`quench`) |
| `first_file_number` | `initial_output` | 0 | Número do primeiro arquivo de saída |

#### Condição Inicial

| Parâmetro | Valores | Padrão | Descrição |
|---|---|---|---|
| `ic` | `random` / `homogeneous` / `ic_file` / `cholesteric` | `random` | Tipo de condição inicial |
| `ic_file` | caminho | — | Arquivo CSV de condição inicial (quando `ic ic_file`) |
| `phi_0`, `theta_0` | float (graus) | 0 | Ângulos inicial (modos `homogeneous`/`cholesteric`) |
| `p0_i` | float | 0 | Passo inicial para condição colestérica |

#### Ancoragem de Superfície (por índice de superfície `n`)

| Parâmetro | Descrição |
|---|---|
| `anchoring_type n <tipo>` | Tipo de ancoragem da superfície `n` |
| `W n <valor>` | Força de ancoragem da superfície `n` |
| `phi_s n <valor>` | Ângulo azimutal de ancoragem da superfície `n` |
| `theta_s n <valor>` | Ângulo polar de ancoragem da superfície `n` |

Tipos disponíveis: `homeotropic`, `fg` (Fournier–Galatola), `rp` (Rapini–Papoular), `strong`.

#### Campo Elétrico

| Parâmetro | Aliases | Descrição |
|---|---|---|
| `elecX`, `elecY`, `elecZ` | `Electric_field_x/y/z`, `EF_x/y/z` | Direção do campo (vetor, não precisa ser unitário) |
| `elecA` | `dielectric_anisotropy`, `Aniso_E` | Anisotropia dielétrica |
| `Ei`, `Ef`, `dE` | `initial_E`, `final_E`, `E_variation` | Intensidade inicial, final e passo do campo (modo `electric`) |

> **`Nk` (neighbour_kind) deve ser sempre `1`.** Vizinhos de segunda ou terceira ordem não são suportados nesta versão — qualquer outro valor faz o programa abortar com uma mensagem de erro explicando o motivo.

---

## Geometrias

**`bulk`** — Rede cúbica sem superfícies; todas as moléculas estão no interior.

**`slab`** — Rede com duas superfícies planas (faces inferior e superior em Z). Simula um filme confinado entre dois substratos.

**`sphere`** — Confinamento esférico inscrito na caixa cúbica. Pontos fora do raio da esfera são marcados como inativos (`pt = 0`).

**`custom`** — Geometria definida por um arquivo CSV externo:

```
boundary_file  minha_geometria.csv
```

O arquivo deve conter colunas `x`, `y`, `z` e `pt`, onde `pt = 0` é sítio desativado, `pt = 1` é interior e `pt ≥ 2` identifica diferentes superfícies.

---

## Modos de Evolução

| Modo | Comportamento |
|---|---|
| `thermal` | Varre a temperatura de `Ti` até `Tf` em passos de `dT`. Em cada temperatura: `MCT` passos de termalização + `MCS` passos de medição. Gera um `director_field_<T*100>.csv` por temperatura. |
| `step` | Executa `fn` blocos independentes de relaxação a temperatura fixa (`Ti`). Útil para estatística/estudo de convergência. |
| `quench` | Como `step`, mas cada bloco aplica dois estágios: `MCT` passos em `Ti` seguidos de passos em `Tf`, simulando um resfriamento brusco. |
| `electric` | Varre a intensidade do campo elétrico de `Ei` até `Ef` em passos de `dE`, a `Ti` fixa. |

Qualquer um dos quatro modos aceita o sufixo `GPU` (ex. `thermalGPU`) para rodar no dispositivo CUDA, caso o binário tenha sido compilado com `make` (não `make CPU`).

---

## Potenciais

**Lebwohl–Lasher (`ll`)** — potencial clássico de pares nemáticos, `U_ij = -A·P₂(cos θ_ij)` mais termos de correção (`B1`, `B2`, `C`).

**Gruhn–Hess / GHRL (`ghrl`)** — potencial elástico contínuo discretizado na rede, baseado nas constantes de Frank (`k11`, `k22`, `k33`). Permite anisotropia elástica real e materiais colestéricos via `p0`.

**Splay-Bend / Pear (`pear`)** — potencial para moléculas com forma de pera, com acoplamento splay-bend e polarização (`P`) reportável no `po.dat`.

---

## Condições de Contorno

| Valor | Descrição |
|---|---|
| `free` | Sítios fora da borda são ignorados no cálculo de energia |
| `periodic` | A rede se repete, simulando um sistema infinito |

Cada eixo pode ter condição independente.

---

## Ancoragem de Superfície

Define como as moléculas próximas a uma superfície se alinham. Cada superfície recebe um índice `n`.

Exemplo — slab com ancoragem homeotrópica nas duas faces:

```
anchoring_type  0  homeotropic
W               0  0.5
anchoring_type  1  homeotropic
W               1  0.5
```

Exemplo — ancoragem planar oblíqua:

```
anchoring_type  0  fg
W               0  0.3
phi_s           0  45.0
theta_s         0  90.0
```

---

## Campo Elétrico

Aplicado como contribuição adicional à energia de cada sítio: `U_elec = -elecA · (n · Ê)²`, onde `Ê` é a direção normalizada `(elecX, elecY, elecZ)`. Para varrer o campo, use `evol electric` com `Ei`, `Ef` e `dE`.

---

## Saídas

Toda simulação **limpa automaticamente** as saídas de execuções anteriores (`ic.csv`, `po.dat`, `director_field_*.csv`) antes de começar — não há risco de resultados antigos se misturarem com os novos.

### `po.dat`

```
T       S        varS      E        varE
0.80    0.712    0.00031   -1.234   0.00018
```

`T` (ou intensidade de campo/índice, conforme o modo) · `S` parâmetro de ordem global · `varS`/`varE` variâncias · `E` energia média por sítio válido. Uma coluna `P` (polarização) é adicionada quando `potential = pear`.

### `ic.csv`
Snapshot da condição inicial, antes do início da evolução.

### `director_field_<id>.csv`

```
x,y,z,nx,ny,nz,S,pt
0,0,0,0.123,-0.456,0.781,0.712,1
```

`(x,y,z)` coordenadas · `(nx,ny,nz)` componentes do diretor · `S` parâmetro de ordem local · `pt` tipo do ponto (`1` interior, `≥2` superfície, `0` inativo).

> Em pontos com `pt = 0` (fora da geometria simulada), `S` é reportado como `NaN` — não há material ali, então não existe parâmetro de ordem a calcular. Ferramentas como pandas, Excel e Origin reconhecem `NaN` automaticamente e o excluem de gráficos e estatísticas.

---

## Estrutura do Projeto

```
McLiCS-backend/
├── include/                       # Cabeçalhos
│   ├── parameters.h               # Struct com todos os parâmetros de simulação
│   ├── evolve.h / evolve.cuh      # Classes de evolução — CPU e GPU
│   ├── geometry.h / geometry.cuh  # Classes de geometria
│   ├── potential.h / potential.cuh
│   ├── anchoring.h                # Tipos de ancoragem de superfície
│   ├── io.h                       # Leitura de parâmetros e escrita de saída
│   ├── ic.h / simulator.h
│   └── parameter_order.h          # Cálculo do tensor de ordem Q
│
├── src/
│   ├── main.cpp                   # Ponto de entrada, banner, limpeza de saídas antigas
│   ├── io.cpp                     # Parser do arquivo de parâmetros
│   ├── simulator.cpp / simulatorGPU.cpp
│   │
│   ├── evolve_thermal.cpp/.cu     # Varredura de temperatura
│   ├── evolve_step.cpp/.cu        # Blocos independentes a T fixa
│   ├── evolve_quench.cpp/.cu      # Resfriamento brusco
│   ├── evolve_electric.cpp/.cu    # Varredura de campo elétrico
│   │
│   ├── geometry_bulk.cpp
│   ├── geometry_slab.cpp
│   ├── geometry_sphere.cpp
│   ├── geometry_custom.cpp
│   │
│   ├── potential.cpp/.cu          # Potenciais LL, GHRL e Pear
│   ├── anchoring_*.cpp            # homeotropic, fg, rp, strong
│   ├── parameter_order.cpp        # Diagonalização do tensor Q (parâmetro de ordem)
│   └── ic.cpp                     # Geração de condição inicial
│
└── Makefile
```

---

## Observações

- O simulador utiliza apenas vizinhos mais próximos (**`Nk = 1`**); vizinhos de segunda ou terceira ordem não são suportados, e o programa aborta com mensagem de erro caso `Nk` diferente de `1` seja especificado.
- A escolha entre CPU e GPU é feita pelo **sufixo do parâmetro `evol`** (ex. `thermal` vs. `thermalGPU`), não por flag de linha de comando — o binário gerado por `make` roda os dois; o gerado por `make CPU` só reconhece os modos sem sufixo `GPU`.
- Em CPU, o OpenMP utiliza todos os núcleos disponíveis, com particionamento estático do trabalho entre threads — o que garante que a mesma simulação, rodada duas vezes na mesma máquina, produza resultados estatisticamente equivalentes.
- O parâmetro de ordem `S` é calculado via diagonalização do tensor de ordem Q de Landau–de Gennes (GSL `gsl_eigen_symmv`).
