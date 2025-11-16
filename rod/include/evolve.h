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
struct dim3 {
  unsigned int x; // Alterado de uint para unsigned int para consistência
  unsigned int y;
  unsigned int z;
};
#endif

class Evolve {
 public:
  Evolve(unsigned int nx, unsigned int ny, unsigned int nz, Geometry *geo, int points)
      : Nx(nx), Ny(ny), Nz(nz), geometry(geo), VallidPoints(points) {}
      
  Evolve() = default; // Construtor padrão opcional

  virtual ~Evolve() = default; 

  virtual int run() = 0; // Tornando a classe base abstrata
  virtual float latice_Potential() const = 0; // Tornando a classe base abstrata e const
  
  unsigned int Nx, Ny, Nz;
  Geometry *geometry;
  int VallidPoints;

  void check_Points(int *pt, Parameters params);
  
  virtual void tester() const {
    printf("here\n");
    fflush(stdout);
  }
};

class EvolveN : public Evolve {
 public:
  EvolveN(float *ni_ptr, int *pt_ptr, Parameters *params_ptr) 
      : Evolve(params_ptr->Nx, params_ptr->Ny, params_ptr->Nz, nullptr, 0), // Inicializa a base
        ni(ni_ptr), pt(pt_ptr), params(params_ptr) {
  };

  int run() override { return 0; };
  float latice_Potential() const override { return 0; };

  void Monte_Carlo_Step(float &ang_var, gsl_rng **r);
  float energy_calculator_GPU() const;
  float energy_calculator() const;
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

  float *ni, *d_ni; 
  int *pt;
  Parameters *params;
  
  // Destrutor virtual para garantir desalocação de ponteiros d_* (device)
  virtual ~EvolveN() {
  }
};

class thermalEvolveN : public EvolveN {
 public:
  thermalEvolveN(float *ni, int *pt, Parameters *params) : EvolveN(ni, pt, params) {};
  int run() override; 
  ~thermalEvolveN() override; 
  
};

class stepEvolveN : public EvolveN {
 public:
  stepEvolveN(float *ni, int *pt, Parameters *params) : EvolveN(ni, pt, params) {};
  int run() override;
  ~stepEvolveN() override;
};

class quenchEvolveN : public EvolveN {
 public:
  quenchEvolveN(float *ni, int *pt, Parameters *params) : EvolveN(ni, pt, params) {};
  int run() override;
  ~quenchEvolveN() override;
};

class electricEvolveN : public EvolveN {
 public:
  electricEvolveN(float *ni, int *pt, Parameters *params) : EvolveN(ni, pt, params) {};
  int run() override;
  ~electricEvolveN() override;
};

#endif