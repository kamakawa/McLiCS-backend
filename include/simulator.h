#ifndef SIMULATOR_H_
#define SIMULATOR_H_

// --- System Includes ---
#include <iostream>
#include <vector>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/evolve_strategy.h"
#include "../include/geometry_strategy.h"
#include "../include/anchoring_strategy.h"
#include "../include/parameters.h"

class simulator {
 public:
  simulator(Parameters *params);
  ~simulator();

  void Setup_simmulation(Parameters &params);
  int run_evolution();

  // Wrapper de snapshot: mantém a API do simulator,
  // mas delega para a implementação oficial em io.cpp
  int print_n(char *fname, Parameters *params);

  // Getters
  GeometryStrategy* get_geometry() { return geometry_strategy; }
  float* get_ni() { return ni; }
  int* get_pt() { return pt; }
  float* get_surface_normals() { return surface_normals; }
  std::vector<AnchoringStrategy*>& get_anchoring_strategies() { return anchoring_strategies; }

 private:
  // Estado do sistema
  float *ni;
  int *pt;
  Parameters *params;

  // Estratégias
  GeometryStrategy *geometry_strategy;
  EvolveStrategy *evolve_strategy;
  std::vector<AnchoringStrategy*> anchoring_strategies;

  // Dados auxiliares
  float *surface_normals;
  int Nx, Ny, Nz;

  // Métodos de configuração
  void setup_evolution_strategy();
  void setup_geometry_strategy(Parameters &params);
  void setup_anchoring_strategies(Parameters &params);
  void setup_potential(Parameters &params);
};

#endif