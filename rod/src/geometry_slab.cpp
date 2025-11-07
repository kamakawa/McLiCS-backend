#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>
#include <memory>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "../include/anchoring.h"
#include "../include/geometry.h"
#include "../include/parameters.h"
#include "../include/potential.h"

Slab_Geometry::Slab_Geometry(int *pt, Parameters *params) : Geometry(params) {
  printf("Geometry: Slab\n");
  nSurfaces = 2;
  
  ns = std::make_unique<float[]>(Nx * Ny * Nz * 3);
  for (int i = 0; i < Nx * Ny * Nz * 3; i++) {
    ns[i] = 0.0f;
  }
  
  surfaces.clear();
  pt = set_point_type_normals(pt, params);

  params->ZBoundtype = "free";
  
  if (params->XBoundtype == "free")
    params->XBound = &Free_Boundary;
  else if (params->XBoundtype == "periodic")
    params->XBound = &Periodic_Boundary;
  else {
    fprintf(stderr, "X boundary condition: %s not implemented \n", params->XBoundtype.c_str());
    exit(2);
  }

  if (params->YBoundtype == "free")
    params->YBound = &Free_Boundary;
  else if (params->YBoundtype == "periodic")
    params->YBound = &Periodic_Boundary;
  else {
    fprintf(stderr, "Y boundary condition: %s not implemented \n", params->YBoundtype.c_str());
    exit(2);
  }

  printf("xbound  %s\n", params->XBoundtype.c_str());
  printf("ybound  %s\n", params->YBoundtype.c_str());
  printf("\n");
}

int *Slab_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  for (int ii = 0; ii < Nx; ii++) {
    for (int jj = 0; jj < Ny; jj++) {
      int kk = 0;
      pt[ii + Nx * (jj + Ny * kk)] = 2;
      ns[(ii + Nx * (jj + Ny * kk)) * 3 + 0] = 0;
      ns[(ii + Nx * (jj + Ny * kk)) * 3 + 1] = 0;
      ns[(ii + Nx * (jj + Ny * kk)) * 3 + 2] = -1;
      
      for (kk = 1; kk < Nz - 1; kk++) {
        pt[ii + Nx * (jj + Ny * kk)] = 0;
      }
      
      kk = Nz - 1;
      pt[ii + Nx * (jj + Ny * kk)] = 3;
      ns[(ii + Nx * (jj + Ny * kk)) * 3 + 0] = 0;
      ns[(ii + Nx * (jj + Ny * kk)) * 3 + 1] = 0;
      ns[(ii + Nx * (jj + Ny * kk)) * 3 + 2] = 1;
    }
  }
  return pt;
}

float Slab_Geometry::lattice_Potential(const nni fullni[7]) {
  float E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};

  E = Geometry::newman_neighbours(fullni);

  if (params->neighbourKind > 1)
    E += Geometry::second_nerghbours(fullni);
  if (params->neighbourKind == 3)
    E += Geometry::third_nerghbours(fullni);
  
  float s[3] = {fullni[7].x, fullni[7].y, fullni[7].z};
  if (fullni[0].pt > 1)
    E += surfaces[fullni[0].pt - 2]->surface_potential(ni, s);
  
  if (params->elecA != 0) 
    E += Potential::Electric_Potential(ni, *params);

  return E;
}