#include <gsl/gsl_rng.h>
#include <string.h>

#include <iostream>

#include "../include/define.h"
//#include "../include/evolve.cuh"
#include "../include/evolve.h"
#include "../include/ic.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"
#include "../include/simulator.h"

simulator::simulator(Parameters *params) 
    // CORREÇÃO: Acesso às dimensões da grade aninhadas em 'lattice'
    : Nx(params->lattice.Nx), Ny(params->lattice.Ny), Nz(params->lattice.Nz) {
  this->params = params;
}

void simulator::Setup_simmulation(Parameters &params) {

  // CORREÇÃO: Acesso às dimensões da grade aninhadas em 'lattice'
  int Nn = 3 * params.lattice.Nx * params.lattice.Ny * params.lattice.Nz;
  ni = (float *)calloc(Nn, sizeof(float));
  pt = (int *)calloc(Nn, sizeof(int));

  if (strcasecmp(params.mc.evol, "thermal") == 0) {
    evolve = new thermalEvolveN(ni, pt, &params);
  } else if (strcasecmp(params.mc.evol, "step") == 0) {
    evolve = new stepEvolveN(ni, pt, &params);
  } else if (strcasecmp(params.mc.evol, "quench") == 0) {
    evolve = new quenchEvolveN(ni, pt, &params);
  } else if (strcasecmp(params.mc.evol, "electric") == 0) {
    evolve = new electricEvolveN(ni, pt, &params);
  } else if (strcasecmp(params.mc.evol, "thermalGPU") == 0) {
  //  evolve = new thermalEvolveNGPU(ni, pt, &params);
  } else if (strcasecmp(params.mc.evol, "stepGPU") == 0) {
  //  evolve = new stepEvolveNGPU(ni, pt, &params);
  } else if (strcasecmp(params.mc.evol, "quenchGPU") == 0) {
  //  evolve = new quenchEvolveNGPU(ni, pt, &params);
  } else if (strcasecmp(params.mc.evol, "electricGPU") == 0) {
  //  evolve = new electricEvolveNGPU(ni, pt, &params);
  } else {

    printf("Evolve %s not defined, try thermal or step\n", params.mc.evol);
    exit(2);
  }

  if (strcasecmp(params.lattice.geometry, "bulk") == 0) {
    evolve->geometry = new Bulk_Geometry(pt, &params);
  } else if (strcasecmp(params.lattice.geometry, "slab") == 0) {
    evolve->geometry = new Slab_Geometry(pt, &params);
  } else if (strcasecmp(params.lattice.geometry, "sphere") == 0) {
    evolve->geometry = new Sphere_Geometry(pt, &params);
  } else if (strcasecmp(params.lattice.geometry, "custom") == 0) {
    evolve->geometry = new Custom_Geometry(pt, &params);
  } else {
    // CORREÇÃO: Acesso a 'geometry' aninhado em 'lattice'
    fprintf(stderr, "Geometry %s not defined\n", params.lattice.geometry);
    exit(2);
  }

  // CORREÇÃO DE ESCOPO: Funções de Boundary estão no escopo global (não em Potential::)
  if (strcasecmp(params.lattice.XBoundtype, "free") == 0)
    params.lattice.XBound = &Free_Boundary; // CORREÇÃO: Removido Potential::
  else if (strcasecmp(params.lattice.XBoundtype, "periodic") == 0)
    params.lattice.XBound = &Periodic_Boundary; // CORREÇÃO: Removido Potential::
  else {
    // CORREÇÃO: Acesso a XBoundtype aninhado em 'lattice'
    fprintf(stderr, "X boundary condition: %s not implemented \n", params.lattice.XBoundtype);
    exit(2);
  }

  if (strcasecmp(params.lattice.YBoundtype, "free") == 0)
    params.lattice.YBound = &Free_Boundary; // CORREÇÃO: Removido Potential::
  else if (strcasecmp(params.lattice.YBoundtype, "periodic") == 0)
    params.lattice.YBound = &Periodic_Boundary; // CORREÇÃO: Removido Potential::
  else {
    // CORREÇÃO: Acesso a YBoundtype aninhado em 'lattice'
    fprintf(stderr, "Y boundary condition: %s not implemented \n", params.lattice.YBoundtype);
    exit(2);
  }

  if (strcasecmp(params.lattice.ZBoundtype, "free") == 0)
    params.lattice.ZBound = &Free_Boundary; // CORREÇÃO: Removido Potential::
  else if (strcasecmp(params.lattice.ZBoundtype, "periodic") == 0)
    params.lattice.ZBound = &Periodic_Boundary; // CORREÇÃO: Removido Potential::
  else {
    // CORREÇÃO: Acesso a ZBoundtype aninhado em 'lattice'
    fprintf(stderr, "Z boundary condition: %s not implemented \n", params.lattice.ZBoundtype);
    exit(2);
  }
  
  evolve->geometry->Boundary_Init(&params);
  evolve->check_Points(pt, params);
  apply_Initial_Condidions(ni, pt, params);


  if (strcasecmp(params.potential.potential, "ll") == 0 || strcasecmp(params.potential.potential, "lebwohl-lahser") == 0) {
    evolve->geometry->bulk_potential = &Potential::Bulk_Energy_Lebwohl_Lasher;
    printf("Using lebwohl-lasher potential\n");
    // CORREÇÃO: Acesso a 'potential' e Correção Lógica de * para ||
  } else if (strcasecmp(params.potential.potential, "ghrl") == 0 || strcasecmp(params.potential.potential, "grun-hess") == 0) {
    evolve->geometry->bulk_potential = &Potential::Bulk_Energy_GHRL;
    // CORREÇÃO: Chamada de função do namespace IO
    IO::setGHRL(params);
    printf("Using gruhn-hess potential\n");
  } else if (strcasecmp(params.potential.potential, "pear") == 0) {
    evolve->geometry->bulk_potential = &Potential::Bulk_Energy_Selinger_Pear;
    printf("Using splay-bend potential\n");
  } else {
    // CORREÇÃO: Acesso a 'potential'
    printf("%s potential not programed.\n Try lebwohl-lasher(LL), pear, BC or gruhn-hess(GHRL)", params.potential.potential);
    exit(2);
  }
}

int simulator::print_n(char *fname, Parameters *params) {
  FILE *output = fopen(fname, "w");
  if (output == 0) {
    perror(fname);
    return 1;
  }
  fprintf(output, "x,y,z,nx,ny,nz,s,pt\n");

  // CORREÇÃO: Acesso às dimensões da grade aninhadas em 'lattice'
  for (int k = 0; k < params->lattice.Nz; k++) {
    for (int j = 0; j < params->lattice.Ny; j++) {
      for (int i = 0; i < params->lattice.Nx; i++) {
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