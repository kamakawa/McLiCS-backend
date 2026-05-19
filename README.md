# McLiCS — Monte Carlo Liquid Crystal Simulator

> Simulador Monte Carlo de cristais líquidos nemáticos desenvolvido como projeto de Iniciação Científica.

**Autor:** Eric Kamakawa  
**Orientadores:** Prof. Rafael Zola e Prof. Rodolfo Teixeira  
**Curso:** Engenharia da Computação – UTFPR, Apucarana

---

## Sumário

1. [Visão Geral](#visão-geral)
2. [Dependências](#dependências)
3. [Compilação](#compilação)
4. [Uso](#uso)
5. [Arquivo de Parâmetros](#arquivo-de-parâmetros)
6. [Geometrias](#geometrias)
7. [Modos de Evolução](#modos-de-evolução)
8. [Potenciais](#potenciais)
9. [Condições de Contorno](#condições-de-contorno)
10. [Ancoragem de Superfície](#ancoragem-de-superfície)
11. [Campo Elétrico](#campo-elétrico)
12. [Saídas](#saídas)
13. [Estrutura do Projeto](#estrutura-do-projeto)

---

## Visão Geral

O McLiCS simula o comportamento de moléculas de cristal líquido nemático em uma rede cúbica tridimensional utilizando o método de Monte Carlo. O simulador suporta execução tanto em **CPU** (via OpenMP) quanto em **GPU** (via CUDA), permitindo ao usuário escolher o dispositivo de acordo com o hardware disponível.

A cada passo Monte Carlo, o simulador calcula a energia de interação entre vizinhos próximos de cada sítio da rede e aceita ou rejeita perturbações angulares das moléculas com base no critério de Metropolis. Ao final de cada temperatura ou passo, são calculados e salvos o parâmetro de ordem escalar **S** e a energia média **E** do sistema.

---

## Dependências

| Dependência | Versão mínima | Necessário para |
|---|---|---|
| GCC | 9+ | Compilação CPU |
| NVCC (CUDA Toolkit) | 11+ | Compilação GPU |
| OpenMP | — | Paralelismo CPU |
| GSL (GNU Scientific Library) | 2.x | Geração de números aleatórios |
| C++ Standard Library | C++17 | Filesystem, chrono |

Instalação das dependências no Ubuntu/Debian:

```bash
sudo apt install g++ libgsl-dev libomp-dev
# Para GPU:
# instale o CUDA Toolkit conforme sua versão de driver em https://developer.nvidia.com/cuda-downloads
```

---

## Compilação

O projeto usa um `Makefile` que gera um único binário unificado (`mc_sim`) capaz de rodar tanto em CPU quanto em GPU.

```bash
make
```

O binário gerado é `./mc_sim`. Internamente ele detecta, em tempo de execução, se deve usar CPU ou GPU com base no arquivo de parâmetros ou nas flags da linha de comando.

---

## Uso

```bash
./mc_sim <arquivo_de_parâmetros> [--gpu | --cpu]
```

**Exemplos:**

```bash
# Usar o dispositivo definido no arquivo de parâmetros
./mc_sim params.txt

# Forçar GPU independente do arquivo de parâmetros
./mc_sim params.txt --gpu

# Forçar CPU independente do arquivo de parâmetros
./mc_sim params.txt --cpu
```

As flags `--gpu` e `--cpu` têm prioridade sobre o campo `device` no arquivo de parâmetros.

---

## Arquivo de Parâmetros

O arquivo de parâmetros é um arquivo de texto com um parâmetro por linha, no formato `chave valor`. Linhas que começam com `#` são comentários e são ignoradas. Os nomes dos parâmetros **não são case-sensitive**.

Exemplo completo:

```
# Rede
Nx          30
Ny          30
Nz          30

# Dispositivo: cpu ou gpu
device      cpu

# Geometria: bulk | slab | sphere | custom
geometry    bulk

# Condições de contorno: free | periodic
xbound      periodic
ybound      periodic
zbound      periodic

# Potencial: ll | ghrl | pear
potential   ll

# Parâmetros do potencial ll
A           1.0

# Modo de evolução: thermal | step | quench | electric
evol        thermal

# Temperatura inicial, final e passo
Ti          1.2
Tf          0.0
dT         -0.02

# Passos Monte Carlo de termalização e medição
MCT         500
MCS         200

# Condição inicial: random | uniform | ic_file
ic          random
```

---

### Referência Completa de Parâmetros

#### Rede e Dispositivo

| Parâmetro | Valores | Padrão | Descrição |
|---|---|---|---|
| `Nx` | inteiro | 16 | Tamanho da rede em X |
| `Ny` | inteiro | 16 | Tamanho da rede em Y |
| `Nz` | inteiro | 16 | Tamanho da rede em Z |
| `device` | `cpu` / `gpu` | `cpu` | Dispositivo de execução |
| `geometry` | `bulk` / `slab` / `sphere` / `custom` | `bulk` | Geometria da simulação |
| `xbound` | `free` / `periodic` | `free` | Condição de contorno em X |
| `ybound` | `free` / `periodic` | `free` | Condição de contorno em Y |
| `zbound` | `free` / `periodic` | `free` | Condição de contorno em Z |

#### Potencial e Constantes Elásticas

| Parâmetro | Valores | Padrão | Descrição |
|---|---|---|---|
| `potential` | `ll` / `ghrl` / `pear` | `ll` | Potencial de interação |
| `A` | float | 1.0 | Parâmetro do potencial LL |
| `B1`, `B2` | float | 0.04, 0.4 | Parâmetros do potencial LL com birrefringência |
| `C` | float | 0.3 | Parâmetro do potencial LL |
| `k11` | float | 1.0 | Constante elástica de splay (GHRL) |
| `k22` | float | 1.0 | Constante elástica de twist (GHRL) |
| `k33` | float | 1.0 | Constante elástica de bend (GHRL) |
| `p0` | float | 0 | Passo do colesterol (GHRL, 0 = nemático) |

#### Evolução e Monte Carlo

| Parâmetro | Aliases | Padrão | Descrição |
|---|---|---|---|
| `evol` | `mode` | `thermal` | Modo de evolução |
| `Ti` | — | 1.2 | Temperatura inicial |
| `Tf` | — | 0.0 | Temperatura final |
| `dT` | — | -0.02 | Passo de temperatura |
| `MCS` | — | 100 | Passos de medição por temperatura |
| `MCT` | — | 100 | Passos de termalização por temperatura |
| `fn` | — | 1 | Número de arquivos (modos step/quench) |
| `first_file_number` | `initial_output` | 0 | Número do primeiro arquivo de saída |

#### Condição Inicial

| Parâmetro | Valores | Padrão | Descrição |
|---|---|---|---|
| `ic` | `random` / `uniform` / `ic_file` | `random` | Tipo de condição inicial |
| `ic_file` | caminho | — | Arquivo CSV de condição inicial (quando `ic ic_file`) |
| `phi_0` | float (graus) | 0 | Ângulo azimutal inicial (modo uniform) |
| `theta_0` | float (graus) | 0 | Ângulo polar inicial (modo uniform) |
| `p0_i` | float | 0 | Passo inicial para condição helicoidal |

#### Ancoragem de Superfície (por índice de superfície `n`)

| Parâmetro | Descrição |
|---|---|
| `anchoring_type n <tipo>` | Tipo de ancoragem da superfície `n` |
| `W n <valor>` | Força de ancoragem da superfície `n` |
| `phi_s n <valor>` | Ângulo azimutal de ancoragem da superfície `n` |
| `theta_s n <valor>` | Ângulo polar de ancoragem da superfície `n` |

Tipos de ancoragem: `homeotropic`, `fg` (Fournier-Galatola), `rp` (Rapini-Papoular), `strong`.

#### Campo Elétrico

| Parâmetro | Aliases | Descrição |
|---|---|---|
| `elecX` | `Electric_field_x`, `EF_x` | Componente X da direção do campo |
| `elecY` | `Electric_field_y`, `EF_y` | Componente Y da direção do campo |
| `elecZ` | `Electric_field_z`, `EF_z` | Componente Z da direção do campo |
| `elecA` | `dielectric_anisotropy`, `Aniso_E` | Anisotropia dielétrica |
| `Ei` | `initial_E` | Intensidade inicial do campo (modo electric) |
| `Ef` | `final_E` | Intensidade final do campo (modo electric) |
| `dE` | `E_variation` | Passo da intensidade do campo |

> **Nota:** `nk` (neighbour_kind) deve ser sempre `1`. Qualquer outro valor causa erro e encerra o programa.

---

## Geometrias

### `bulk`
Rede cúbica sem superfícies. Todas as moléculas estão no interior.

### `slab`
Rede com duas superfícies planas (faces Z inferior e superior). Usada para simular filmes confinados entre dois substratos.

### `sphere`
Rede com confinamento esférico. Moléculas fora da esfera são tratadas como pontos de fronteira.

### `custom`
Geometria definida por um arquivo CSV externo. O arquivo é especificado com:

```
boundary_file  minha_geometria.csv
```

O arquivo CSV deve conter colunas `x`, `y`, `z` e `pt`, onde `pt = 0` indica sítio desativado, `pt = 1` é interior e `pt >= 2` identifica diferentes superfícies.

---

## Modos de Evolução

### `thermal`
Varre a temperatura de `Ti` até `Tf` com passo `dT`. Em cada temperatura, realiza `MCT` passos de termalização seguidos de `MCS` passos de medição. Gera um arquivo `director_field_<T*100>.csv` por temperatura.

### `step`
Realiza `fn` blocos independentes de relaxação a temperatura fixa (`Ti`). Cada bloco gera um arquivo de saída numerado a partir de `first_file_number`. Útil para coletar estatísticas ou estudar convergência.

### `quench`
Similar ao `step`, mas cada bloco aplica dois estágios: `MCT` passos a `Ti` (alta temperatura) seguidos de `MCT * dT` passos a `Tf` (baixa temperatura). Simula um resfriamento brusco (*quench*).

### `electric`
Varre a intensidade do campo elétrico de `Ei` até `Ef` com passo `dE`, a temperatura fixa `Ti`. Requer os parâmetros de campo elétrico configurados.

---

## Potenciais

### Lebwohl-Lasher (`ll`)
Potencial clássico de interação entre pares de moléculas nemáticas:

```
U_ij = -A * P2(cos θ_ij)
```

onde `P2` é o polinômio de Legendre de segunda ordem e `θ_ij` é o ângulo entre os diretores `n_i` e `n_j`. Parâmetros: `A`, `B1`, `B2`, `C`.

### Gruhn-Hess (GHRL) (`ghrl`)
Potencial elástico contínuo discretizado na rede, baseado nas constantes de Frank (`k11`, `k22`, `k33`). Permite simular nemáticos com anisotropia elástica real e materiais colestéricos (via `p0`).

### Splay-Bend / Pear (`pear`)
Potencial para moléculas com forma de pera (*pear-shaped*), que exibem acoplamento splay-bend.

---

## Condições de Contorno

| Valor | Descrição |
|---|---|
| `free` | Condição de superfície livre — sítios fora da borda são ignorados no cálculo de energia |
| `periodic` | Condição periódica — a rede se repete, simulando um sistema infinito |

Cada eixo (X, Y, Z) pode ter condição independente.

---

## Ancoragem de Superfície

A ancoragem define como as moléculas próximas a uma superfície se alinham. Cada superfície recebe um índice `n` (começando em 0 para geometria `slab`, ou conforme definido no arquivo de boundary para `custom`).

Exemplo para geometria slab com ancoragem homeotrófica nas duas faces:

```
anchoring_type  0  homeotropic
W               0  0.5
anchoring_type  1  homeotropic
W               1  0.5
```

Exemplo para ancoragem planar oblíqua:

```
anchoring_type  0  fg
W               0  0.3
phi_s           0  45.0
theta_s         0  90.0
```

---

## Campo Elétrico

O campo elétrico é aplicado como uma contribuição adicional à energia de cada sítio, proporcional à anisotropia dielétrica `elecA` e à intensidade `elecE`:

```
U_elec = -elecA * (n · E)^2
```

A direção do campo é definida pelo vetor unitário `(elecX, elecY, elecZ)`.

Para varrer o campo, use o modo `electric` e configure `Ei`, `Ef` e `dE`.

---

## Saídas

### `po.dat`
Arquivo de texto com os observáveis medidos a cada passo principal (temperatura, intensidade de campo ou índice de passo). Formato:

```
T      S      varS      E      varE
0.80   0.712  0.00031   -1.234  0.00018
...
```

- `T` — temperatura (ou intensidade de campo/índice, dependendo do modo)
- `S` — parâmetro de ordem escalar médio
- `varS` — variância de S
- `E` — energia média por sítio
- `varE` — variância de E

### `ic.csv`
Snapshot da condição inicial, gerado antes do início da evolução.

### `director_field_<id>.csv`
Snapshot do campo diretor ao fim de cada passo principal. Formato:

```
x,y,z,nx,ny,nz,S,pt
0,0,0,0.123,-0.456,0.781,0.712,1
...
```

- `(x, y, z)` — coordenadas do sítio
- `(nx, ny, nz)` — componentes do diretor
- `S` — parâmetro de ordem local
- `pt` — tipo do ponto (`1` = interior, `≥2` = superfície, `0` = inativo)

> Arquivos de saída de execuções anteriores são **removidos automaticamente** no início de cada simulação.

---

## Estrutura do Projeto

```
backend_refactored/
├── include/                  # Cabeçalhos
│   ├── parameters.h          # Struct com todos os parâmetros de simulação
│   ├── progress.h            # Barra de progresso, temporizador e ETA
│   ├── evolve.h / evolve.cuh # Classes de evolução CPU e GPU
│   ├── geometry.h / .cuh     # Classes de geometria
│   ├── potential.h / .cuh    # Funções de energia
│   ├── anchoring.h           # Ancoragem de superfície
│   ├── io.h                  # Leitura de parâmetros e escrita de saída
│   └── ...
│
├── src/                      # Implementações
│   ├── main.cpp              # Ponto de entrada, banner, roteamento CPU/GPU
│   ├── io.cpp                # Parser do arquivo de parâmetros
│   │
│   ├── evolve_thermal.cpp/.cu    # Evolução por varredura de temperatura
│   ├── evolve_step.cpp/.cu       # Evolução por passos fixos
│   ├── evolve_quench.cpp/.cu     # Evolução com quench térmico
│   ├── evolve_electric.cpp/.cu   # Evolução por varredura de campo elétrico
│   │
│   ├── geometry_bulk.cpp     # Geometria bulk (sem confinamento)
│   ├── geometry_slab.cpp     # Geometria de filme plano
│   ├── geometry_sphere.cpp   # Geometria esférica
│   ├── geometry_custom.cpp   # Geometria definida por arquivo externo
│   │
│   ├── potential.cpp/.cu     # Potenciais LL, GHRL e Pear
│   ├── anchoring_*.cpp       # Tipos de ancoragem de superfície
│   ├── ic.cpp                # Geração de condição inicial
│   └── ...
│
└── Makefile
```

---

## Observações

- O simulador utiliza vizinhos mais próximos (`nk = 1`). Vizinhos de segunda ou terceira ordem não são suportados — o programa aborta com mensagem de erro se `nk` diferente de `1` for especificado no arquivo de parâmetros.
- Em modo GPU, o número de threads é determinado automaticamente pela arquitetura da placa. Em modo CPU, o OpenMP utiliza todos os núcleos disponíveis.
- O parâmetro de ordem `S` é calculado via diagonalização do tensor de ordem Q de Landau-de Gennes.
