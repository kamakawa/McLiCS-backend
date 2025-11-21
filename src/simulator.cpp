#include "../include/simulator.h"

// --- System Includes ---
#include <gsl/gsl_rng.h>
#include <string.h>
#include <iostream>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/evolve_strategy.h"  
#include "../include/ic.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

// Construtor
simulator::simulator(Parameters *params) 
    : Nx(params->Nx), Ny(params->Ny), Nz(params->Nz) {
  this->params = params;
  
  // Inicializa ponteiros como nulos para seguranca
  this->ni = nullptr;
  this->pt = nullptr;
  this->geometry = nullptr;       
  this->evolve_strategy = nullptr; 
}

// Destrutor (Limpeza de Memoria)
simulator::~simulator() {
  if (this->evolve_strategy != nullptr) {
    delete this->evolve_strategy;
  }

  if (this->geometry != nullptr) {
    delete this->geometry;
  }

  if (this->ni != nullptr) free(this->ni);
  if (this->pt != nullptr) free(this->pt);
}

void simulator::Setup_simmulation(Parameters &params) {
  int Nn = 3 * params.Nx * params.Ny * params.Nz;
  ni = (float *)calloc(Nn, sizeof(float));
  pt = (int *)calloc(Nn, sizeof(int));

  if (strcasecmp(params.geometry, "bulk") == 0) {
    geometry = new Bulk_Geometry(pt, &params);
  } else if (strcasecmp(params.geometry, "slab") == 0) {
    geometry = new Slab_Geometry(pt, &params);
  } else if (strcasecmp(params.geometry, "sphere") == 0) {
    geometry = new Sphere_Geometry(pt, &params);
  } else if (strcasecmp(params.geometry, "custom") == 0) {
    geometry = new Custom_Geometry(pt, &params);
  } else {
    fprintf(stderr, "Geometry %s not defined\n", params.geometry);
    exit(2);
  }

  // --- Condicao de Contorno X ---
  if (strcasecmp(params.XBoundtype, "free") == 0)
    params.XBound = &Potential::Free_Boundary;
  else if (strcasecmp(params.XBoundtype, "periodic") == 0)
    params.XBound = &Potential::Periodic_Boundary;
  else {
    fprintf(stderr, "X boundary condition: %s not implemented \n", params.XBoundtype);
    exit(2);
  }

  // --- Condicao de Contorno Y ---
  if (strcasecmp(params.YBoundtype, "free") == 0)
    params.YBound = &Potential::Free_Boundary;
  else if (strcasecmp(params.YBoundtype, "periodic") == 0)
    params.YBound = &Potential::Periodic_Boundary;
  else {
    fprintf(stderr, "Y boundary condition: %s not implemented \n", params.YBoundtype);
    exit(2);
  }

  // --- Condicao de Contorno Z ---
  if (strcasecmp(params.ZBoundtype, "free") == 0)
    params.ZBound = &Potential::Free_Boundary;
  else if (strcasecmp(params.ZBoundtype, "periodic") == 0)
    params.ZBound = &Potential::Periodic_Boundary;
  else {
    fprintf(stderr, "Z boundary condition: %s not implemented \n", params.ZBoundtype);
    exit(2);
  }

  // Inicializacao de Geometria e Pontos
  geometry->Boundary_Init(&params);
  
  // Aplicacao das Condicoes Iniciais
  apply_Initial_Condidions(ni, pt, params);

  // --- Selecao do Potencial de Interacao ---
  if ((strcasecmp(params.potential, "ll") == 0) || (strcasecmp(params.potential, "lebwohl-lahser") == 0)) {
    geometry->bulk_potential = &Potential::Bulk_Energy_Lebwohl_Lasher;
    printf("Using lebwohl-lasher potential\n");
  } 
  else if ((strcasecmp(params.potential, "ghrl") == 0) || (strcasecmp(params.potential, "grun-hess") == 0)) {
    geometry->bulk_potential = &Potential::Bulk_Energy_GHRL;
    setGHRL(params);
    printf("Using gruhn-hess potential\n");
  } 
  else if (strcasecmp(params.potential, "pear") == 0) {
    geometry->bulk_potential = &Potential::Bulk_Energy_Selinger_Pear;
    printf("Using splay-bend potential\n");
  } 
  else {
    printf("%s potential not programed.\n Try lebwohl-lasher(LL), pear, BC or gruhn-hess(GHRL)", params.potential);
    exit(2);
  }

  setup_evolution_strategy();
}

void simulator::setup_evolution_strategy() {
  evolve_strategy = EvolveStrategyFactory::create(params);
  
  if (!evolve_strategy) {
    std::cerr << "Error: Could not create evolution strategy!" << std::endl;
    exit(1);
  }
  
  printf("Evolution strategy configured: %s\n", 
         evolve_strategy->getName().c_str());
}

int simulator::run_evolution() {
  if (!evolve_strategy) {
    std::cerr << "Error: No evolution strategy configured!" << std::endl;
    return -1;
  }
  
  if (!geometry) {
    std::cerr << "Error: No geometry configured!" << std::endl;
    return -1;
  }
  
  printf("\n=== Starting Evolution ===\n");
  printf("Strategy: %s\n", evolve_strategy->getName().c_str());
  printf("System size: %dx%dx%d\n", params->Nx, params->Ny, params->Nz);
  printf("Temperature: Ti=%g, Tf=%g\n", params->Ti, params->Tf);
  printf("Monte Carlo: MCT=%d, MCS=%d\n\n", params->MCT, params->MCS);
  
  // Executa a estratégia
  return evolve_strategy->run(ni, pt, params, geometry);
}

// Imprime snapshot (Wrapper para a funcao io)
int simulator::print_n(char *fname, Parameters *params) {
  FILE *output = fopen(fname, "w");
  if (output == 0) {
    perror(fname);
    return 1;
  }
  fprintf(output, "x,y,z,nx,ny,nz,s,pt\n");

  for (int k = 0; k < params->Nz; k++) {
    for (int j = 0; j < params->Ny; j++) {
      for (int i = 0; i < params->Nx; i++) {
        fprintf(output, "%d,%d,%d,%g,%g,%g,1,%d\n", i, j, k,
                nix(i, j, k), niy(i, j, k), niz(i, j, k), pti(i, j, k));
      }
    }
  }
  printf("Snapshot saved in %s\n", fname);
  fflush(stdout);
  fclose(output);
  return 0;
}