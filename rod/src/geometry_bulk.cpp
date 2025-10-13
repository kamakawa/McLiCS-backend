#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "../include/potential.h" 
#include "../include/anchoring.h"
#include "../include/geometry.h"
#include "../include/parameters.h"

// Note: A linha 'using namespace Potential;' não é necessária se você usar 'Potential::',
// mas garante que o compilador encontre as funções. No seu caso, vamos garantir que
// todas as referências sigam a nova arquitetura de forma consistente.

Bulk_Geometry::Bulk_Geometry(int *pt, Parameters *params) : Geometry(params) {
  printf("Geometry: Bulk\n");
  // CORREÇÃO 1, 2, 3: Acesso às dimensões da grade aninhadas em 'lattice'
  ns = (float *)calloc(params->lattice.Nx * params->lattice.Ny * params->lattice.Nz * 3, sizeof(float));

  nSurfaces = 0;
  surfaces = std::vector<class Anchoring *>(nSurfaces);
  pt = set_point_type_normals(pt, params);

  // CORREÇÃO 4, 5: Acesso a XBoundtype e XBound aninhados em 'lattice'
  if (strcasecmp(params->lattice.XBoundtype, "free") == 0)
    params->lattice.XBound = Potential::Free_Boundary; // CORREÇÃO: Namespace Potential
  else if (strcasecmp(params->lattice.XBoundtype, "periodic") == 0)
    params->lattice.XBound = Potential::Periodic_Boundary; // CORREÇÃO: Namespace Potential
  else {
    fprintf(stderr, "X boundary condition: %s not implemented \n", params->lattice.XBoundtype);
    exit(2);
  }

  // CORREÇÃO 6, 7: Acesso a YBoundtype e YBound aninhados em 'lattice'
  if (strcasecmp(params->lattice.YBoundtype, "free") == 0)
    params->lattice.YBound = &Potential::Free_Boundary; // CORREÇÃO: Namespace Potential
  else if (strcasecmp(params->lattice.YBoundtype, "periodic") == 0)
    params->lattice.YBound = &Potential::Periodic_Boundary; // CORREÇÃO: Namespace Potential
  else {
    fprintf(stderr, "Y boundary condition: %s not implemented \n", params->lattice.YBoundtype);
    exit(2);
  }

  // CORREÇÃO 8, 9: Acesso a ZBoundtype e ZBound aninhados em 'lattice'
  if (strcasecmp(params->lattice.ZBoundtype, "free") == 0)
    params->lattice.ZBound = &Potential::Free_Boundary; // CORREÇÃO: Namespace Potential
  else if (strcasecmp(params->lattice.ZBoundtype, "periodic") == 0)
    params->lattice.ZBound = &Potential::Periodic_Boundary; // CORREÇÃO: Namespace Potential
  else {
    fprintf(stderr, "Z boundary condition: %s not implemented \n", params->lattice.ZBoundtype);
    exit(2);
  }

  // CORREÇÃO 10, 11, 12: Acesso aos Boundarytypes aninhados em 'lattice'
  printf("xbound  %s\n", params->lattice.XBoundtype);
  printf("ybound  %s\n", params->lattice.YBoundtype);
  printf("zbound  %s\n", params->lattice.ZBoundtype);
  printf("\n");
}

int *Bulk_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  // CORREÇÃO 13: Acesso às dimensões da grade aninhadas em 'lattice'
  for (int ii = 0; ii < params->lattice.Nx * params->lattice.Ny * params->lattice.Nz; ii++) pt[ii] = 1;

  return pt;
}

float Bulk_Geometry::latice_Potential(const nni fullni[7]) {
  float rij[3];
  float E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
  E = Geometry::newman_neighbours(fullni);
  
  // CORREÇÃO 14: Acesso a neighbourKind aninhado em 'neighbourhood'
  if (params->neighbourhood.neighbourKind > 1)
    E += Geometry::second_nerghbours(fullni);
    
  // CORREÇÃO 15: Acesso a neighbourKind aninhado em 'neighbourhood'
  if (params->neighbourhood.neighbourKind == 3)
    E += Geometry::third_nerghbours(fullni);
    
  // CORREÇÃO 16, 17: Acesso a elecA aninhado em 'electric' e Namespace Potential
  if (params->electric.elecA!=0) 
    E+=Potential::Electric_Potential(ni,params);

  return E;
}