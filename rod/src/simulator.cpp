#include <gsl/gsl_rng.h>
#include <string.h>
#include <memory>
#include <algorithm>
#include <cstdio>
#include <string>

#include <iostream>

#include "../include/define.h"
#include "../include/evolve.h"
#include "../include/ic.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"
#include "../include/simulator.h"

extern int Free_Boundary(int &ii, int NN);
extern int Periodic_Boundary(int &ii, int NN);

namespace PotentialWrapper {
    float Bulk_Energy_Lebwohl_Lasher_Wrapper(float ni[3], float nj[3], Parameters *params, float rij[3], int nk) {
        return Potential::Bulk_Energy_Lebwohl_Lasher(ni, nj, *params, rij, nk);
    }
    float Bulk_Energy_GHRL_Wrapper(float ni[3], float nj[3], Parameters *params, float rij[3], int nk) {
        return Potential::Bulk_Energy_GHRL(ni, nj, *params, rij, nk);
    }
    float Bulk_Energy_Selinger_Pear_Wrapper(float ni[3], float nj[3], Parameters *params, float rij[3], int nk) {
        return Potential::Bulk_Energy_Selinger_Pear(ni, nj, *params, rij, nk);
    }
}

simulator::simulator(Parameters *params) : Nx(params->Nx), Ny(params->Ny), Nz(params->Nz) {
  this->params = params;
}

void simulator::Setup_simmulation(Parameters &params) {
  int Nn = params.Nx * params.Ny * params.Nz;
  int Nn3 = 3 * Nn;
  
  ni = std::make_unique<float[]>(Nn3); 
  pt = std::make_unique<int[]>(Nn);    
  
  std::fill(ni.get(), ni.get() + Nn3, 0.0f);
  std::fill(pt.get(), pt.get() + Nn, 0);
  
  float *ni_raw = ni.get();
  int *pt_raw = pt.get();

  if (params.evol == "thermal") {
    evolve = std::make_unique<thermalEvolveN>(ni_raw, pt_raw, &params, this);
  } else if (params.evol == "step") {
    evolve = std::make_unique<stepEvolveN>(ni_raw, pt_raw, &params, this);
  } else if (params.evol == "quench") {
    evolve = std::make_unique<quenchEvolveN>(ni_raw, pt_raw, &params, this);
  } else if (params.evol == "electric") {
    evolve = std::make_unique<electricEvolveN>(ni_raw, pt_raw, &params, this);
  } else {
    printf("Evolve %s not defined, try thermal or step\n", params.evol.c_str());
    exit(2);
  }

  if (params.geometry == "bulk") {
    evolve->geometry = new Bulk_Geometry(pt_raw, &params);
  } else if (params.geometry == "slab") {
    evolve->geometry = new Slab_Geometry(pt_raw, &params);
  } else if (params.geometry == "sphere") {
    evolve->geometry = new Sphere_Geometry(pt_raw, &params);
  } else if (params.geometry == "custom") {
    evolve->geometry = new Custom_Geometry(pt_raw, &params);
  } else {
    fprintf(stderr, "Geometry %s not defined\n", params.geometry.c_str());
    exit(2);
  }

  if (params.XBoundtype == "free")
    params.XBound = &Free_Boundary;
  else if (params.XBoundtype == "periodic")
    params.XBound = &Periodic_Boundary;
  else {
    fprintf(stderr, "X boundary condition: %s not implemented \n", params.XBoundtype.c_str());
    exit(2);
  }

  if (params.YBoundtype == "free")
    params.YBound = &Free_Boundary;
  else if (params.YBoundtype == "periodic")
    params.YBound = &Periodic_Boundary;
  else {
    fprintf(stderr, "Y boundary condition: %s not implemented \n", params.YBoundtype.c_str());
    exit(2);
  }

  if (params.ZBoundtype == "free")
    params.ZBound = &Free_Boundary;
  else if (params.ZBoundtype == "periodic")
    params.ZBound = &Periodic_Boundary;
  else {
    fprintf(stderr, "Z boundary condition: %s not implemented \n", params.ZBoundtype.c_str());
    exit(2);
  }

  evolve->geometry->Boundary_Init(&params);
  evolve->check_Points(pt_raw, params);
  
  apply_Initial_Condidions(ni, pt, params);

  if (params.potential == "ll" || params.potential == "lebwohl-lahser") {
    evolve->geometry->bulk_potential = &PotentialWrapper::Bulk_Energy_Lebwohl_Lasher_Wrapper;
    printf("Using lebwohl-lasher potential\n");
  } else if (params.potential == "ghrl" || params.potential == "grun-hess") {
    evolve->geometry->bulk_potential = &PotentialWrapper::Bulk_Energy_GHRL_Wrapper;
    IO::setGHRL(params);
    printf("Using gruhn-hess potential\n");
  } else if (params.potential == "pear") {
    evolve->geometry->bulk_potential = &PotentialWrapper::Bulk_Energy_Selinger_Pear_Wrapper;
    printf("Using splay-bend potential\n");
  } else {
    printf("%s potential not programed.\n Try lebwohl-lasher(LL), pear, BC or gruhn-hess(GHRL)", params.potential.c_str());
    exit(2);
  }
}

int simulator::print_n(const std::string& fname, Parameters &params) {
  return IO::print_n(fname, ni, params, pt);
}