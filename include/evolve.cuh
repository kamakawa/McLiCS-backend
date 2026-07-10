#ifndef EVOL_CUH_
#define EVOL_CUH_
#include "../include/potential.h"
#include "../include/geometry.h"
#include <iostream>
#include "../include/parameters.h"
#include "../include/define.h"
#include "../include/evolve.h"

class EvolveNGPU: public Evolve{
public:

EvolveNGPU(float *ni, int *pt, Parameters *params):ni(ni),pt(pt),params(params),Nx(params->Nx),Ny(params->Ny),Nz(params->Nz){};
virtual int run(){return 0;};
float energy_calculator_GPU();
void Monte_Carlo_Step_GPU( float &ang_var, gsl_rng * r);
__device__ void MC_GPU(dim3 tick, int di, int dj, int dk, float *ni, int *pt, float *d_T, float ang_var, curandState *const rngStates, unsigned int *acceptance, Parameters *params);
protected:
dim3 tick;
curandState *d_rngStates = 0;
int *d_pt=0;
float *d_T;
unsigned int *d_acc=0;
Parameters *d_params=0;
int Nx, Ny, Nz;
float *ni, *d_ni;
int *pt;
Parameters *params;

};
class thermalEvolveNGPU: public EvolveNGPU{
  public:
    thermalEvolveNGPU(float *ni, int *pt, Parameters *params);
    int run();
    ~thermalEvolveNGPU();
  int Nx, Ny, Nz;
  float *ni;
  int *pt;
};

class stepEvolveNGPU: public EvolveNGPU{
public:
  stepEvolveNGPU(float *ni, int *pt, Parameters *params);
  int run();
  ~stepEvolveNGPU();
int Nx, Ny, Nz;
float *ni;
int *pt;
};
class quenchEvolveNGPU: public EvolveNGPU{
  public:
    quenchEvolveNGPU(float *ni, int *pt, Parameters *params);
    int run();
    ~quenchEvolveNGPU();
  int Nx, Ny, Nz;
  float *ni;
  int *pt;
};
class electricEvolveNGPU: public EvolveNGPU{
  public:
    electricEvolveNGPU(float *ni, int *pt, Parameters *params);
    int run();
    ~electricEvolveNGPU();
  int Nx, Ny, Nz;
  float *ni;
  int *pt;
};
__global__ void initRNG(curandState *const rngStates, const unsigned int seed, dim3 tick);


#endif
