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

Sphere_Geometry::Sphere_Geometry(int *pt, Parameters *params) : Geometry(params) {
  printf("Geometry: Sphere\n");
  nSurfaces = 1;
  
  params->XBoundtype = "free";
  params->YBoundtype = "free";
  params->ZBoundtype = "free";
  
  ns = std::make_unique<float[]>(Nx * Ny * Nz * 3);
  for (int i = 0; i < Nx * Ny * Nz * 3; i++) {
    ns[i] = 0.0f;
  }
  
  surfaces.clear();
  pt = set_point_type_normals(pt, params);
  printf("\n");
}

int *Sphere_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  float Hx = Nx / 2.0f;
  float Hy = Ny / 2.0f;
  float Hz = Nz / 2.0f;
  
  for (int ii = 0; ii < Nx; ii++) {
    for (int jj = 0; jj < Ny; jj++) {
      for (int kk = 0; kk < Nz; kk++) {
        float Rx = ii - Hx;
        float Ry = jj - Hy;
        float Rz = kk - Hz;
        float radius = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);
        
        if (radius < Hz - 1)
          pt[ii + Nx * (jj + Ny * kk)] = 0;
        else if (radius < Hz + 1) {
          pt[ii + Nx * (jj + Ny * kk)] = 2;
          ns[(ii + Nx * (jj + Ny * kk)) * 3 + 0] = Rx / radius;
          ns[(ii + Nx * (jj + Ny * kk)) * 3 + 1] = Ry / radius;
          ns[(ii + Nx * (jj + Ny * kk)) * 3 + 2] = Rz / radius;
        } else
          pt[ii + Nx * (jj + Ny * kk)] = -1;
      }
    }
  }
  return pt;
}

float Sphere_Geometry::lattice_Potential(const nni fullni[7]) {
  float E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
  float s[3] = {fullni[7].x, fullni[7].y, fullni[7].z};
  
  E = Geometry::newman_neighbours(fullni);

  if (params->neighbourKind > 1)
    E += Geometry::second_nerghbours(fullni);
  if (params->neighbourKind == 3)
    E += Geometry::third_nerghbours(fullni);
  
  if (fullni[0].pt > 1)
    E += surfaces[fullni[0].pt - 2]->surface_potential(ni, s);
  
  if (params->elecA != 0) 
    E += Potential::Electric_Potential(ni, *params);

  return E;
}