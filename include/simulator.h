#ifndef SIMULATOR_H_
#define SIMULATOR_H_

// --- System Includes ---
#include <iostream>
#include <vector>  // ← NOVO

// --- Project Includes ---
#include "../include/define.h"
#include "../include/evolve_strategy.h"
#include "../include/geometry_strategy.h"    // ← NOVO
#include "../include/anchoring_strategy.h"  // ← NOVO
#include "../include/parameters.h"

class simulator {
 public:
  simulator(Parameters *params);
  ~simulator();

  void Setup_simmulation(Parameters &params);
  int run_evolution();
  int print_n(char *fname, Parameters *params);
  
  // Getters atualizados
  GeometryStrategy* get_geometry() { return geometry_strategy; }  // ← NOVO
  float* get_ni() { return ni; }
  int* get_pt() { return pt; }
  float* get_surface_normals() { return surface_normals; }  // ← NOVO
  std::vector<AnchoringStrategy*>& get_anchoring_strategies() { return anchoring_strategies; }  // ← NOVO

 private:
  // Estado do sistema
  float *ni;
  int *pt;
  Parameters *params;
  
  // Estratégias
  GeometryStrategy *geometry_strategy;           // ← NOVO (substitui Geometry*)
  EvolveStrategy *evolve_strategy;
  std::vector<AnchoringStrategy*> anchoring_strategies;  // ← NOVO
  
  // Dados auxiliares
  float *surface_normals;  // ← NOVO
  int Nx, Ny, Nz;
  
  // Métodos de configuração
  void setup_evolution_strategy();
  void setup_geometry_strategy(Parameters &params);      // ← NOVO
  void setup_anchoring_strategies(Parameters &params);   // ← NOVO
  void setup_potential(Parameters &params);              // ← NOVO
};

#endif