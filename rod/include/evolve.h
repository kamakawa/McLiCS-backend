//Modifiquei
//Removidos membros duplicados Nx, Ny, Nz, ni, pt das classes derivadas
//Adicionado destrutor virtual em Evolve
//Destrutores nas classes derivadas com override = default
#ifndef EVOL_H_
#define EVOL_H_
#include <iostream>

#include "../include/define.h"
#include "../include/geometry.h"
#include "../include/parameters.h"
#include "../include/potential.h"
#ifdef CUDA__
// #pragma message ( "Cuda Compilation" )
#include <cuda_runtime.h>
#include <curand_kernel.h>
#else
// #pragma message ( "CPU Compilation" )
struct dim3 {
  uint x;
  uint y;
  uint z;
};
#endif

class Evolve {
 public:
  // NOVO: Destrutor virtual para garantir que a desalocação de memória seja correta (sem memory leaks).
  virtual ~Evolve() = default; 
  
  // As funções virtuais não precisam ser alteradas na base, mas se fossem puras, estariam aqui.
  virtual int run() { return 0; }; 
  virtual float latice_Potential() { return 0; };
  
  unsigned int Nx, Ny, Nz;
  Geometry *geometry;
  int VallidPoints;

  void check_Points(int *pt, Parameters params);
  virtual void tester() {
    printf("here\n");
    fflush(stdout);
  }
};

class EvolveN : public Evolve {
 public:
  // Construtor com a correção da leitura de 'Parameters' (acesso à struct 'lattice')
  EvolveN(float *ni, int *pt, Parameters *params)
   : ni(ni), pt(pt), params(params), Nx(params->lattice.Nx), Ny(params->lattice.Ny), Nz(params->lattice.Nz){};

  // Implementação da função run() com a palavra-chave 'override'
  virtual int run() override { return 0; }; 
  
  void Monte_Carlo_Step(float &ang_var, gsl_rng **r);
  float energy_calculator_GPU();
  float energy_calculator();
  void Monte_Carlo_Step_GPU(float &ang_var, gsl_rng *r);

 protected:
  dim3 tick;
#ifdef CUDA__
  curandState *d_rngStates = 0;
#endif
  int *d_pt = 0;
  float *d_T;
  unsigned int *d_acc = 0;
  Parameters *d_params = 0;
  
  // Estas variáveis já são herdadas de Evolve (ou deveriam estar lá).
  // Se elas forem mantidas aqui, o compilador deve estar ciente de que as dimensões
  // Nx, Ny, Nz que você inicializou acima pertencem à EvolveN.
  // Vamos remover a declaração delas nas classes derivadas:
  int Nx, Ny, Nz;
  float *ni, *d_ni;
  int *pt;
  Parameters *params;
};

class thermalEvolveN : public EvolveN {
 public:
  thermalEvolveN(float *ni, int *pt, Parameters *params);
  int run() override; // Uso do 'override'
  // CORREÇÃO: Usando `= default` para o destrutor (conforme seu comentário), e removendo membros duplicados.
  ~thermalEvolveN() = default; 
  // Removendo: int Nx, Ny, Nz; float *ni; int *pt;
};

class stepEvolveN : public EvolveN {
 public:
  stepEvolveN(float *ni, int *pt, Parameters *params);
  int run() override; // Uso do 'override'
  // CORREÇÃO: Usando `= default` para o destrutor, e removendo membros duplicados.
  ~stepEvolveN() = default; 
  // Removendo: int Nx, Ny, Nz; float *ni; int *pt;
};

class quenchEvolveN : public EvolveN {
 public:
  quenchEvolveN(float *ni, int *pt, Parameters *params);
  int run() override; // Uso do 'override'
  // CORREÇÃO: Usando `= default` para o destrutor, e removendo membros duplicados.
  ~quenchEvolveN() = default; 
  // Removendo: int Nx, Ny, Nz; float *ni; int *pt;
};

class electricEvolveN : public EvolveN {
 public:
  electricEvolveN(float *ni, int *pt, Parameters *params);
  int run() override; // Uso do 'override'
  // CORREÇÃO: Usando `= default` para o destrutor, e removendo membros duplicados.
  ~electricEvolveN() = default; 
  // Removendo: int Nx, Ny, Nz; float *ni; int *pt;
};

#endif