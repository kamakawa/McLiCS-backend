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
  virtual int run() = 0;  //Método puro
  virtual float lattice_Potential() = 0;  //Método puro
  virtual ~Evolve() = default;  //Destrutor virtual

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
  EvolveN(float *ni, int *pt, Parameters *params) : ni(ni), pt(pt), params(params) {
    //Usar membros da classe base
    this->Nx = params->Nx;
    this->Ny = params->Ny;
    this->Nz = params->Nz;
  };
  
  //Destrutor virtual para limpeza GPU
  virtual ~EvolveN() {
    cleanup_gpu_memory();
  }
  
  virtual int run() override { return 0; };
  virtual float lattice_Potential() override { return 0; };
  
  void Monte_Carlo_Step(float &ang_var, gsl_rng **r);
  float energy_calculator_GPU();
  float energy_calculator();
  void Monte_Carlo_Step_GPU(float &ang_var, gsl_rng *r);

 protected:
  
  void cleanup_gpu_memory() {
#ifdef CUDA__
    if(d_rngStates) {
      cudaFree(d_rngStates);
      d_rngStates = 0;
    }
    if(d_pt) {
      cudaFree(d_pt);
      d_pt = 0;
    }
    if(d_ni) {
      cudaFree(d_ni);
      d_ni = 0;
    }
    if(d_params) {
      cudaFree(d_params);
      d_params = 0;
    }
    if(d_acc) {
      cudaFree(d_acc);
      d_acc = 0;
    }
    if(d_T) {
      cudaFree(d_T);
      d_T = 0;
    }
#endif
  }

  dim3 tick;
#ifdef CUDA__
  curandState *d_rngStates = 0;
#endif
  int *d_pt = 0;
  float *d_T = 0;
  unsigned int *d_acc = 0;
  Parameters *d_params = 0;
  
  float *ni, *d_ni = 0;
  int *pt;
  Parameters *params;
};

class thermalEvolveN : public EvolveN {
 public:
  thermalEvolveN(float *ni, int *pt, Parameters *params) : EvolveN(ni, pt, params) {};
  int run() override;
  ~thermalEvolveN() override = default;
};

class stepEvolveN : public EvolveN {
 public:
  stepEvolveN(float *ni, int *pt, Parameters *params) : EvolveN(ni, pt, params) {};
  int run() override;
  ~stepEvolveN() override = default;
};

class quenchEvolveN : public EvolveN {
 public:
  quenchEvolveN(float *ni, int *pt, Parameters *params) : EvolveN(ni, pt, params) {};
  int run() override;
  ~quenchEvolveN() override = default;
};

class electricEvolveN : public EvolveN {
 public:
  electricEvolveN(float *ni, int *pt, Parameters *params) : EvolveN(ni, pt, params) {};
  int run() override;
  ~electricEvolveN() override = default;
};

#endif