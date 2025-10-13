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

Slab_Geometry::Slab_Geometry(int *pt, Parameters *params) : Geometry(params) {
  printf("Geometry: Slab\n");
  nSurfaces = 2;
  // Alteração 1, 2, 3: Acesso às dimensões da grade aninhadas em 'lattice'
  ns = (float *)calloc(params->lattice.Nx * params->lattice.Ny * params->lattice.Nz * 3, sizeof(float));
  surfaces = std::vector<class Anchoring *>(nSurfaces);
  pt = set_point_type_normals(pt, params);

  // Alteração 4: Acesso a ZBoundtype aninhado em 'lattice'
  sprintf(params->lattice.ZBoundtype, "free");
  
  // Alteração 5, 6: Acesso a XBoundtype e XBound aninhados em 'lattice'
  if (strcasecmp(params->lattice.XBoundtype, "free") == 0)
    params->lattice.XBound = &Potential::Free_Boundary; // Alteração 7: Chamada de função do namespace Potential
  else if (strcasecmp(params->lattice.XBoundtype, "periodic") == 0)
    params->lattice.XBound = &Potential::Periodic_Boundary; // Alteração 8: Chamada de função do namespace Potential
  else {
    fprintf(stderr, "X boundary condition: %s not implemented \n", params->lattice.XBoundtype);
    exit(2);
  }

  // Alteração 9, 10: Acesso a YBoundtype e YBound aninhados em 'lattice'
  if (strcasecmp(params->lattice.YBoundtype, "free") == 0)
    params->lattice.YBound = &Potential::Free_Boundary; // Alteração 11: Chamada de função do namespace Potential
  else if (strcasecmp(params->lattice.YBoundtype, "periodic") == 0)
    params->lattice.YBound = &Potential::Periodic_Boundary; // Alteração 12: Chamada de função do namespace Potential
  else {
    fprintf(stderr, "Y boundary condition: %s not implemented \n", params->lattice.YBoundtype);
    exit(2);
  }

  // Alteração 13, 14: Acesso a XBoundtype e YBoundtype aninhados em 'lattice'
  printf("xbound  %s\n", params->lattice.XBoundtype);
  printf("ybound  %s\n", params->lattice.YBoundtype);
  printf("\n");
}

int *Slab_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  int kk;
  // Alteração 15, 16, 17: Acesso às dimensões da grade aninhadas em 'lattice'
  for (int ii = 0; ii < params->lattice.Nx; ii++) {
    for (int jj = 0; jj < params->lattice.Ny; jj++) {
      kk = 0;
      pt[ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)] = 2;
      ns[(ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)) * 3 + 0] = 0;
      ns[(ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)) * 3 + 1] = 0;
      ns[(ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)) * 3 + 2] = -1;
      // Alteração 18, 19: Acesso a Nz aninhado em 'lattice'
      for (kk = 1; kk < params->lattice.Nz - 1; kk++) {
        pt[ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)] = 1;
      }
      // Alteração 20, 21: Acesso a Nz aninhado em 'lattice'
      kk = params->lattice.Nz - 1;
      pt[ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)] = 3;
      ns[(ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)) * 3 + 0] = 0;
      ns[(ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)) * 3 + 1] = 0;
      ns[(ii + params->lattice.Nx * (jj + params->lattice.Ny * kk)) * 3 + 2] = 1;
    }
  }
  return pt;
}

float Slab_Geometry::latice_Potential(const nni fullni[7]) {
  float rij[3];
  float E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};

  E = Geometry::newman_neighbours(fullni);

  // Alteração 22: Acesso a neighbourKind aninhado em 'neighbourhood'
  if (params->neighbourhood.neighbourKind > 1)
    E += Geometry::second_nerghbours(fullni);
    
  // Alteração 23: Acesso a neighbourKind aninhado em 'neighbourhood'
  if (params->neighbourhood.neighbourKind == 3)
    E += Geometry::third_nerghbours(fullni);
    
  float s[3] = {fullni[7].x, fullni[7].y, fullni[7].z};
  if (fullni[0].pt > 1)
    E += surfaces[fullni[0].pt - 2]->surface_potential(ni, s);
    
  // Alteração 24: Acesso a elecA aninhado em 'electric'
  // Alteração 25: Chamada de função do namespace Potential
  if (params->electric.elecA!=0) 
    E+=Potential::Electric_Potential(ni,params);

  return E;
}