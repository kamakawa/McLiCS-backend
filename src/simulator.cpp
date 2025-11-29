#include "../include/simulator.h"

// --- System Includes ---
#include <gsl/gsl_rng.h>
#include <string.h>
#include <iostream>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/evolve_strategy.h"  
#include "../include/geometry_strategy.h"
#include "../include/anchoring_strategy.h"
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
  this->geometry_strategy = nullptr;       
  this->evolve_strategy = nullptr; 
  this->surface_normals = nullptr;
}

// Destrutor (Limpeza de Memoria)
simulator::~simulator() {
  // Libera estratégia de evolução
  if (this->evolve_strategy != nullptr) {
    delete this->evolve_strategy;
  }

  // Libera estratégia de geometria
  if (this->geometry_strategy != nullptr) {
    delete this->geometry_strategy;
  }

  // Libera estratégias de ancoramento
  for (auto* anchoring : anchoring_strategies) {
    delete anchoring;
  }
  anchoring_strategies.clear();

  // Libera arrays de dados
  if (this->ni != nullptr) free(this->ni);
  if (this->pt != nullptr) free(this->pt);
  if (this->surface_normals != nullptr) free(this->surface_normals);
}

void simulator::Setup_simmulation(Parameters &params) {
  int Nn = 3 * params.Nx * params.Ny * params.Nz;
  ni = (float *)calloc(Nn, sizeof(float));
  pt = (int *)calloc(Nn, sizeof(int));

  // Configura estratégia de geometria
  setup_geometry_strategy(params);

  // Configura condições de contorno através da estratégia
  geometry_strategy->initializeBoundaries(&params);

  // Configura tipos de pontos e normais
  surface_normals = (float *)calloc(params.Nx * params.Ny * params.Nz * 3, sizeof(float));
  pt = geometry_strategy->setPointTypes(pt, &params, surface_normals);

  // Configura estratégias de ancoramento
  setup_anchoring_strategies(params);

  // Aplicacao das Condicoes Iniciais
  apply_Initial_Condidions(ni, pt, params);

  // CORREÇÃO CRÍTICA: Configura potencial de interação
  setup_potential(params);

  // Configura estratégia de evolução
  setup_evolution_strategy();
}

// Configura estratégia de geometria
void simulator::setup_geometry_strategy(Parameters &params) {
  geometry_strategy = GeometryStrategyFactory::create(&params);
  
  if (!geometry_strategy) {
    std::cerr << "Error: Could not create geometry strategy for type: " 
              << params.geometry << std::endl;
    exit(1);
  }
  
  printf("Geometry strategy configured: %s\n", 
         geometry_strategy->getName().c_str());
}

// Configura estratégias de ancoramento
void simulator::setup_anchoring_strategies(Parameters &params) {
  int num_surfaces = geometry_strategy->getNumSurfaces();
  
  if (num_surfaces > 0) {
    // Pré-aloca o vector com a capacidade necessária
    anchoring_strategies.resize(num_surfaces);
    
    // Inicializa as estratégias de ancoramento
    AnchoringStrategyFactory::initializeAnchoringStrategies(&params, anchoring_strategies);
    
    printf("Configured %d anchoring strategies\n", num_surfaces);
  }
}

// CORREÇÃO: Configura potencial de interação corretamente
void simulator::setup_potential(Parameters &params) {
    float (*potential_func)(float[3], float[3], Parameters*, float[3], int) = nullptr;
    
    if ((strcasecmp(params.potential, "ll") == 0) || (strcasecmp(params.potential, "lebwohl-lahser") == 0)) {
        potential_func = &Potential::Bulk_Energy_Lebwohl_Lasher;
        printf("Using lebwohl-lasher potential\n");
    } 
    else if ((strcasecmp(params.potential, "ghrl") == 0) || (strcasecmp(params.potential, "grun-hess") == 0)) {
        potential_func = &Potential::Bulk_Energy_GHRL;
        setGHRL(params);
        printf("Using gruhn-hess potential\n");
    } 
    else if (strcasecmp(params.potential, "pear") == 0) {
        potential_func = &Potential::Bulk_Energy_Selinger_Pear;
        printf("Using splay-bend potential\n");
    } 
    else {
        printf("%s potential not programed.\n Try lebwohl-lasher(LL), pear, BC or gruhn-hess(GHRL)", params.potential);
        exit(2);
    }
    
    // CORREÇÃO: Configurar o potencial na estratégia de geometria
    if (geometry_strategy && potential_func) {
        geometry_strategy->setBulkPotential(potential_func);
        printf("Bulk potential configured in geometry strategy\n");
    }
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
  
  if (!geometry_strategy) {
    std::cerr << "Error: No geometry strategy configured!" << std::endl;
    return -1;
  }
  
  printf("\n=== Starting Evolution ===\n");
  printf("Strategy: %s\n", evolve_strategy->getName().c_str());
  printf("Geometry: %s\n", geometry_strategy->getName().c_str());
  printf("System size: %dx%dx%d\n", params->Nx, params->Ny, params->Nz);
  printf("Temperature: Ti=%g, Tf=%g\n", params->Ti, params->Tf);
  printf("Monte Carlo: MCT=%d, MCS=%d\n\n", params->MCT, params->MCS);
  
  // Executa a estratégia passando os novos parâmetros
  return evolve_strategy->run(ni, pt, params, geometry_strategy, 
                             surface_normals, anchoring_strategies);
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