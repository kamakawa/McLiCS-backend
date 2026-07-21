#include <gsl/gsl_rng.h>
#include <string.h>

#include <iostream>

#include "../include/define.h"
#include "../include/evolve.h"
//~ #include "../include/evolve.cuh"
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
}

void simulator::Setup_simmulation(Parameters &params) {
  int Nn = 3 * params.Nx * params.Ny * params.Nz;
  ni = (float *)calloc(Nn, sizeof(float));
  pt = (int *)calloc(Nn, sizeof(int));

  if (strcasecmp(params.evol, "thermal") == 0) {
    evolve = new thermalEvolveN(ni, pt, &params);
  } else if (strcasecmp(params.evol, "step") == 0) {
    evolve = new stepEvolveN(ni, pt, &params);
  } else if (strcasecmp(params.evol, "quench") == 0) {
    evolve = new quenchEvolveN(ni, pt, &params);
  } else if (strcasecmp(params.evol, "electric") == 0) {
    evolve = new electricEvolveN(ni, pt, &params);
  } else {
    fprintf(stderr, "  [ERROR] Evol mode '%s' not implemented. Options: thermal | step | quench | electric\n", params.evol);
    exit(2);
  }

  if (strcasecmp(params.geometry, "bulk") == 0) {
    evolve->geometry = new Bulk_Geometry(pt, &params);
  } else if (strcasecmp(params.geometry, "slab") == 0) {
    evolve->geometry = new Slab_Geometry(pt, &params);
  } else if (strcasecmp(params.geometry, "sphere") == 0) {
    evolve->geometry = new Sphere_Geometry(pt, &params);
  } else if (strcasecmp(params.geometry, "custom") == 0) {
    evolve->geometry = new Custom_Geometry(pt, &params);
  } else {
    fprintf(stderr, "  [ERROR] Geometry '%s' not implemented. Options: bulk | slab | sphere | custom\n", params.geometry);
    exit(2);
  }

  if (strcasecmp(params.XBoundtype, "free") == 0)
    params.XBound = &Free_Boundary;
  else if (strcasecmp(params.XBoundtype, "periodic") == 0)
    params.XBound = &Periodic_Boundary;
  else {
    fprintf(stderr, "  [ERROR] X boundary '%s' not implemented. Options: free | periodic\n", params.XBoundtype);
    exit(2);
  }

  if (strcasecmp(params.YBoundtype, "free") == 0)
    params.YBound = &Free_Boundary;
  else if (strcasecmp(params.YBoundtype, "periodic") == 0)
    params.YBound = &Periodic_Boundary;
  else {
    fprintf(stderr, "  [ERROR] Y boundary '%s' not implemented. Options: free | periodic\n", params.YBoundtype);
    exit(2);
  }

  if (strcasecmp(params.ZBoundtype, "free") == 0)
    params.ZBound = &Free_Boundary;
  else if (strcasecmp(params.ZBoundtype, "periodic") == 0)
    params.ZBound = &Periodic_Boundary;
  else {
    fprintf(stderr, "  [ERROR] Z boundary '%s' not implemented. Options: free | periodic\n", params.ZBoundtype);
    exit(2);
  }

  evolve->geometry->Boundary_Init(&params);
  evolve->check_Points(pt, params);
  apply_Initial_Condidions(ni, pt, params);

  if (strcasecmp(params.potential, "ll") * strcasecmp(params.potential, "lebwohl-lahser") == 0) {
    evolve->geometry->bulk_potential = &Bulk_Energy_Lebwohl_Lasher;
    print_field("Potential:", "Lebwohl-Lasher");
    printf("\n");
  } else if (strcasecmp(params.potential, "ghrl") * strcasecmp(params.potential, "grun-hess") == 0) {
    evolve->geometry->bulk_potential = &Bulk_Energy_GHRL;
    setGHRL(params);
    print_field("Potential:", "Gruhn-Hess (GHRL)");
    printf("\n");
  } else if (strcasecmp(params.potential, "pear") == 0) {
    evolve->geometry->bulk_potential = &Bulk_Energy_Selinger_Pear;
    print_field("Potential:", "Splay-bend (Pear)");
    printf("\n");
  } else {
    fprintf(stderr, "  [ERROR] Potential '%s' not implemented. Options: ll | ghrl | pear\n", params.potential);
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

  for (int k = 0; k < params->Nz; k++) {
    for (int j = 0; j < params->Ny; j++) {
      for (int i = 0; i < params->Nx; i++) {
        fprintf(output, "%d,%d,%d,%g,%g,%g,1,%d\n", i, j, k,
                nix(i, j, k), niy(i, j, k), niz(i, j, k), pti(i, j, k));
      }
    }
  }
  printf("  Snapshot -> %s\n", fname);
  fflush(stdout);
  fclose(output);
  return 0;
}
