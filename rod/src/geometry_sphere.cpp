#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <memory> 

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
  
  ns = std::unique_ptr<float[]>(new float[Nx * Ny * Nz * 3]()); 
  
  surfaces.resize(nSurfaces); 
  
  params->XBound = &Free_Boundary;
  params->YBound = &Free_Boundary;
  params->ZBound = &Free_Boundary;
  
  pt = set_point_type_normals(pt, params);
  printf("\n");
}

int *Sphere_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  int ii, jj, kk;
  float Rx, Ry, Rz;
  float Hx = Nx / 2.0f; 
  float Hy = Ny / 2.0f;
  float Hz = Nz / 2.0f;
  
  float *ns_ptr = ns.get(); 
  
  for (ii = 0; ii < Nx; ii++) {
    for (jj = 0; jj < Ny; jj++) {
      for (kk = 0; kk < Nz; kk++) {
        Rx = (float)ii - Hx;
        Ry = (float)jj - Hy;
        Rz = (float)kk - Hz;
        float radius = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);
        if (radius < Hz - 1)
          pt[ii + Nx * (jj + Ny * kk)] = 1;
        else if (radius < Hz + 1) {
          pt[ii + Nx * (jj + Ny * kk)] = 2;
          ns_ptr[(ii + Nx * (jj + Ny * kk)) * 3 + 0] = Rx / radius;
          ns_ptr[(ii + Nx * (jj + Ny * kk)) * 3 + 1] = Ry / radius;
          ns_ptr[(ii + Nx * (jj + Ny * kk)) * 3 + 2] = Rz / radius;
        } else
          pt[ii + Nx * (jj + Ny * kk)] = 0;
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
    
  if (params->elecA!=0) 
    E += Potential::Electric_Potential(ni, *params);

  return E;
}