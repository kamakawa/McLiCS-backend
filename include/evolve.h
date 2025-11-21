#ifndef EVOL_H_
#define EVOL_H_

#include <iostream>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/geometry.h"
#include "../include/parameters.h"
#include "../include/potential.h"

// --- CUDA / CPU Compatibility ---
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

// Classe base abstrata para os processos de evolucao
class Evolve {
 public:
    virtual ~Evolve() = default; // Destrutor virtual para seguranca de memoria

    virtual int run() { return 0; }
    virtual float latice_Potential() { return 0; }

    unsigned int Nx, Ny, Nz;
    Geometry *geometry;
    int VallidPoints;

    void check_Points(int *pt, Parameters params);

    virtual void tester() {
        printf("here\n");
        fflush(stdout);
    }
};

// Classe base intermediaria que gerencia estado comum e logica de Monte Carlo
class EvolveN : public Evolve {
 public:
    EvolveN(float *ni, int *pt, Parameters *params) 
        : ni(ni), pt(pt), params(params), 
          Nx(params->Nx), Ny(params->Ny), Nz(params->Nz) {}

    virtual int run() override { return 0; }

    // Metodos de calculo
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

// Evolucao termica do sistema
class thermalEvolveN : public EvolveN {
 public:
    thermalEvolveN(float *ni, int *pt, Parameters *params);
    ~thermalEvolveN();
    
    int run() override;
    
    int Nx, Ny, Nz;
    float *ni;
    int *pt;
};

// Evolucao por passos discretos
class stepEvolveN : public EvolveN {
 public:
    stepEvolveN(float *ni, int *pt, Parameters *params);
    ~stepEvolveN();

    int run() override;
    
    int Nx, Ny, Nz;
    float *ni;
    int *pt;
};

// Evolucao tipo Quench (resfriamento rapido)
class quenchEvolveN : public EvolveN {
 public:
    quenchEvolveN(float *ni, int *pt, Parameters *params);
    ~quenchEvolveN();

    int run() override;
    
    int Nx, Ny, Nz;
    float *ni;
    int *pt;
};

// Evolucao sob influencia de campo eletrico
class electricEvolveN : public EvolveN {
 public:
    electricEvolveN(float *ni, int *pt, Parameters *params);
    ~electricEvolveN();

    int run() override;
    
    int Nx, Ny, Nz;
    float *ni;
    int *pt;
};

#endif