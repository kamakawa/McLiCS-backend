#include <gsl/gsl_rng.h>
#include <string.h>
#include <iostream>

#include "../include/define.h"
#include "../include/evolve.cuh"
#include "../include/evolve.h"
#include "../include/ic.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"
#include "../include/simulator.h"

simulator::simulator(Parameters *params) : Nx(params->Nx), Ny(params->Ny), Nz(params->Nz) {
  this->params = params;
}

simulator::~simulator() {
  // CORREÇÃO: Liberar memória alocada no construtor
  if (ni != nullptr) {
    free(ni);
    ni = nullptr;
  }
  if (pt != nullptr) {
    free(pt);
    pt = nullptr;
  }
  // CORREÇÃO: Liberar objeto evolve se existir
  if (evolve != nullptr) {
    delete evolve;
    evolve = nullptr;
  }
}

void simulator::Setup_simmulation(Parameters &params) {
  int Nn = 3 * params.Nx * params.Ny * params.Nz;
  ni = (float *)calloc(Nn, sizeof(float));
  pt = (int *)calloc(Nn, sizeof(int));

  // Configuração do tipo de evolução (CPU e GPU)
  if (strcasecmp(params.evol, "thermal") == 0) {
    evolve = new thermalEvolveN(ni, pt, &params);
  } else if (strcasecmp(params.evol, "step") == 0) {
    evolve = new stepEvolveN(ni, pt, &params);
  } else if (strcasecmp(params.evol, "quench") == 0) {
    evolve = new quenchEvolveN(ni, pt, &params);
  } else if (strcasecmp(params.evol, "electric") == 0) {
    evolve = new electricEvolveN(ni, pt, &params);
  } else if (strcasecmp(params.evol, "thermalGPU") == 0) {
    evolve = new thermalEvolveNGPU(ni, pt, &params);
  } else if (strcasecmp(params.evol, "stepGPU") == 0) {
    evolve = new stepEvolveNGPU(ni, pt, &params);
  } else if (strcasecmp(params.evol, "quenchGPU") == 0) {
    evolve = new quenchEvolveNGPU(ni, pt, &params);
  } else if (strcasecmp(params.evol, "electricGPU") == 0) {
    evolve = new electricEvolveNGPU(ni, pt, &params);
  } else {
    printf("Evolve %s not defined, try thermal or step\n", params.evol);
    exit(2);
  }

  // Configuração da geometria
  if (strcasecmp(params.geometry, "bulk") == 0) {
    evolve->geometry = new Bulk_Geometry(pt, &params);
  } else if (strcasecmp(params.geometry, "slab") == 0) {
    evolve->geometry = new Slab_Geometry(pt, &params);
  } else if (strcasecmp(params.geometry, "sphere") == 0) {
    evolve->geometry = new Sphere_Geometry(pt, &params);
  } else if (strcasecmp(params.geometry, "custom") == 0) {
    evolve->geometry = new Custom_Geometry(pt, &params);
  } else {
    fprintf(stderr, "Geometry %s not defined\n", params.geometry);
    exit(2);
  }

  // Configuração das condições de contorno X
  if (strcasecmp(params.XBoundtype, "free") == 0) {
    params.XBound = &Free_Boundary;
  } else if (strcasecmp(params.XBoundtype, "periodic") == 0) {
    params.XBound = &Periodic_Boundary;
  } else {
    fprintf(stderr, "X boundary condition: %s not implemented \n", params.XBoundtype);
    exit(2);
  }

  // Configuração das condições de contorno Y
  if (strcasecmp(params.YBoundtype, "free") == 0) {
    params.YBound = &Free_Boundary;
  } else if (strcasecmp(params.YBoundtype, "periodic") == 0) {
    params.YBound = &Periodic_Boundary;
  } else {
    fprintf(stderr, "Y boundary condition: %s not implemented \n", params.YBoundtype);
    exit(2);
  }

  // Configuração das condições de contorno Z
  if (strcasecmp(params.ZBoundtype, "free") == 0) {
    params.ZBound = &Free_Boundary;
  } else if (strcasecmp(params.ZBoundtype, "periodic") == 0) {
    params.ZBound = &Periodic_Boundary;
  } else {
    fprintf(stderr, "Z boundary condition: %s not implemented \n", params.ZBoundtype);
    exit(2);
  }

  // Inicialização final
  evolve->geometry->Boundary_Init(&params);
  evolve->check_Points(pt, params);
  apply_Initial_Condidions(ni, pt, params);

  // Configuração do potencial
  if (strcasecmp(params.potential, "ll") * strcasecmp(params.potential, "lebwohl-lahser") == 0) {
    evolve->geometry->bulk_potential = &Bulk_Energy_Lebwohl_Lasher;
    printf("Using lebwohl-lasher potential\n");
  } else if (strcasecmp(params.potential, "ghrl") * strcasecmp(params.potential, "grun-hess") == 0) {
    evolve->geometry->bulk_potential = &Bulk_Energy_GHRL;
    setGHRL(params);
    printf("Using gruhn-hess potential\n");
  } else if (strcasecmp(params.potential, "pear") == 0) {
    evolve->geometry->bulk_potential = &Bulk_Energy_Selinger_Pear;
    printf("Using splay-bend potential\n");
  } else {
    printf("%s potential not programed.\n Try lebwohl-lasher(LL), pear, BC or gruhn-hess(GHRL)", params.potential);
    exit(2);
  }
}

int simulator::print_n(char *fname, Parameters *params) {
  // CORREÇÃO: Reutiliza a função já existente e testada do namespace FileIO
  return FileIO::print_n(fname, ni, *params, pt);
}