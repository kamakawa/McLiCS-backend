#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>

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
  // Alteração 1, 2, 3: Acesso aos Boundarytypes aninhados em 'lattice'
  sprintf(params->lattice.XBoundtype, "free");
  sprintf(params->lattice.YBoundtype, "free");
  sprintf(params->lattice.ZBoundtype, "free");
  
  // Alteração 4, 5, 6: Acesso às dimensões da grade aninhadas em 'lattice'
  ns = (float *)calloc(params->lattice.Nx * params->lattice.Ny * params->lattice.Nz * 3, sizeof(float));
  surfaces = std::vector<class Anchoring *>(nSurfaces);
  pt = set_point_type_normals(pt, params);
  printf("\n");
}

int *Sphere_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  int ii, jj, kk;
  float Rx, Ry, Rz;
  
  // Alteração 7, 8, 9: Acesso às dimensões da grade aninhadas em 'lattice'
  float Hx = params->lattice.Nx / 2;
  float Hy = params->lattice.Ny / 2;
  float Hz = params->lattice.Nz / 2;
  
  // Alteração 10, 11, 12: Acesso às dimensões da grade aninhadas em 'lattice'
  for (ii = 0; ii < params->lattice.Nx; ii++) {
    for (jj = 0; jj < params->lattice.Ny; jj++) {
      for (kk = 0; kk < params->lattice.Nz; kk++) {
        Rx = ii - Hx;
        Ry = jj - Hy;
        Rz = kk - Hz;
        float radius = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);
        
        // Alteração 13, 14, 15: Acesso às dimensões da grade aninhadas em 'lattice'
        if (radius < params->lattice.Nz / 2 - 1)
          pt[ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)] = 1;
        else if (radius < params->lattice.Nz / 2 + 1) {
          pt[ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)] = 2;
          ns[(ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)) * 3 + 0] = Rx / radius;
          ns[(ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)) * 3 + 1] = Ry / radius;
          ns[(ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)) * 3 + 2] = Rz / radius;
        } else
          pt[ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)] = 0;
      }
    }
  }
  return pt;
}

float Sphere_Geometry::latice_Potential(const nni fullni[7]) {
  float rij[3];
  float E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
  float s[3] = {fullni[7].x, fullni[7].y, fullni[7].z};
  E = Geometry::newman_neighbours(fullni);

  // Alteração 16: Acesso a neighbourKind aninhado em 'neighbourhood'
  if (params->neighbourhood.neighbourKind > 1)
    E += Geometry::second_nerghbours(fullni);
    
  // Alteração 17: Acesso a neighbourKind aninhado em 'neighbourhood'
  if (params->neighbourhood.neighbourKind == 3)
    E += Geometry::third_nerghbours(fullni);
    
  if (fullni[0].pt > 1)
    E += surfaces[fullni[0].pt - 2]->surface_potential(ni, s);
    
  // Alteração 18: Acesso a elecA aninhado em 'electric'
  // Alteração 19: Chamada de função do namespace Potential
  if (params->electric.elecA!=0) 
    E+=Potential::Electric_Potential(ni,params);

  return E;
}