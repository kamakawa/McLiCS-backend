#ifndef EVOL_H_
#define EVOL_H_
#include <iostream>
#include <memory>

#include "../include/define.h"
#include "../include/geometry.h"
#include "../include/parameters.h"
#include "../include/potential.h"

// Forward declarations para separar concerns CUDA
#ifdef CUDA__
#include <cuda_runtime.h>
#include <curand_kernel.h>
#else
struct dim3 {
  uint x;
  uint y;
  uint z;
};
#endif

class Evolve {
 public:
  virtual ~Evolve() = default;
  virtual int run() = 0;
  virtual float latice_Potential() = 0;
  
  // Getters para acesso controlado
  unsigned int getNx() const { return Nx; }
  unsigned int getNy() const { return Ny; }
  unsigned int getNz() const { return Nz; }
  Geometry* getGeometry() const { return geometry; }
  int getValidPoints() const { return ValidPoints; }
  
  void check_Points(int *pt, Parameters params);
  virtual void tester() {
    printf("here\n");
    fflush(stdout);
  }

 protected:
  Evolve(unsigned int nx, unsigned int ny, unsigned int nz, Geometry* geom) 
    : Nx(nx), Ny(ny), Nz(nz), geometry(geom) {}
    
  unsigned int Nx, Ny, Nz;
  Geometry *geometry;
  int ValidPoints;
};

class EvolveN : public Evolve {
 public:
  EvolveN(float *ni, int *pt, Parameters *params);
  virtual ~EvolveN();
  
  // Interface pública comum
  void Monte_Carlo_Step(float &ang_var, gsl_rng **r);
  float energy_calculator_GPU();
  float energy_calculator();
  void Monte_Carlo_Step_GPU(float &ang_var, gsl_rng *r);
  
  // Getters
  float* getNi() const { return ni; }
  int* getPt() const { return pt; }
  Parameters* getParams() const { return params; }

 protected:
  // Recursos geridos pela classe base
  float *ni, *d_ni;
  int *pt, *d_pt;
  Parameters *params;
  
  // Estado de execução
  dim3 tick;
  
#ifdef CUDA__
  // Recursos CUDA
  curandState *d_rngStates = nullptr;
  float *d_T = nullptr;
  unsigned int *d_acc = nullptr;
  Parameters *d_params = nullptr;
#endif
};

// Classes derivadas simplificadas - sem duplicação de atributos
class thermalEvolveN : public EvolveN {
 public:
  thermalEvolveN(float *ni, int *pt, Parameters *params);
  int run() override;
};

class stepEvolveN : public EvolveN {
 public:
  stepEvolveN(float *ni, int *pt, Parameters *params);
  int run() override;
};

class quenchEvolveN : public EvolveN {
 public:
  quenchEvolveN(float *ni, int *pt, Parameters *params);
  int run() override;
};

class electricEvolveN : public EvolveN {
 public:
  electricEvolveN(float *ni, int *pt, Parameters *params);
  int run() override;
};

#endif