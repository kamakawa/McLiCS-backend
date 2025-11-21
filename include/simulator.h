#ifndef SIMULATOR_H_
#define SIMULATOR_H_

// --- System Includes ---
#include <iostream>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/evolve_strategy.h"  // ← NOVO
#include "../include/geometry.h"
#include "../include/parameters.h"

class simulator {
 public:
  // Construtor e Destrutor
  simulator(Parameters *params);
  ~simulator();

  // Configuracao da simulacao
  void Setup_simmulation(Parameters &params);  
  
  // Execucao da evolucao
  int run_evolution();  
  
  // Metodos de acesso e utilidade
  int print_n(char *fname, Parameters *params);
  Geometry* get_geometry() { return geometry; }
  float* get_ni() { return ni; }
  int* get_pt() { return pt; }

  int Nx, Ny, Nz;
  Geometry *geometry;
  float *ni;
  int *pt;

 private:
  // Estado do sistema
  Parameters *params;
  
  EvolveStrategy *evolve_strategy;
  
  
  // Configuracao da estrategia
  void setup_evolution_strategy();  
};

#endif