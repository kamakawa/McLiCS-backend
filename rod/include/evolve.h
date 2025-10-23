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
  //Destrutor virtual para garantir que a desalocação de memória seja correta
  virtual ~Evolve() = default; 
  
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
  EvolveN(float *ni, int *pt, Parameters *params)
   : ni(ni), pt(pt), params(params), Nx(params->lattice.Nx), Ny(params->lattice.Ny), Nz(params->lattice.Nz){};

  // override indica que este método substitui um método virtual da classe base
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
  
  int Nx, Ny, Nz;
  float *ni, *d_ni;
  int *pt;
  Parameters *params;
};

class thermalEvolveN : public EvolveN {
 public:
  thermalEvolveN(float *ni, int *pt, Parameters *params);
  int run() override; //override
  // destrutor com override
  ~thermalEvolveN() = default; 
};

class stepEvolveN : public EvolveN {
 public:
  stepEvolveN(float *ni, int *pt, Parameters *params);
  int run() override; // override
  ~stepEvolveN() = default; 
};

class quenchEvolveN : public EvolveN {
 public:
  quenchEvolveN(float *ni, int *pt, Parameters *params);
  int run() override; //override
  ~quenchEvolveN() = default; 
};

class electricEvolveN : public EvolveN {
 public:
  electricEvolveN(float *ni, int *pt, Parameters *params);
  int run() override; // override
  ~electricEvolveN() = default; 
};

#endif