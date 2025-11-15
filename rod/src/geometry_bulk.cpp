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

Bulk_Geometry::Bulk_Geometry(int *pt, Parameters *params) : Geometry(params) {
  printf("Geometry: Bulk\n");
  ns = std::make_unique<float[]>(Nx * Ny * Nz * 3);  // CORREÇÃO: smart pointer

  nSurfaces = 0;
  surfaces = std::vector<std::unique_ptr<Anchoring>>(nSurfaces);  // CORREÇÃO: smart pointers
  pt = set_point_type_normals(pt, params);

  // Configuração das condições de contorno X
  if (strcasecmp(params->XBoundtype, "free") == 0) {
    params->XBound = &Free_Boundary;
  } else if (strcasecmp(params->XBoundtype, "periodic") == 0) {
    params->XBound = &Periodic_Boundary;
  } else {
    fprintf(stderr, "X boundary condition: %s not implemented \n", params->XBoundtype);
    exit(2);
  }

  // Configuração das condições de contorno Y
  if (strcasecmp(params->YBoundtype, "free") == 0) {
    params->YBound = &Free_Boundary;
  } else if (strcasecmp(params->YBoundtype, "periodic") == 0) {
    params->YBound = &Periodic_Boundary;
  } else {
    fprintf(stderr, "Y boundary condition: %s not implemented \n", params->YBoundtype);
    exit(2);
  }

  // Configuração das condições de contorno Z
  if (strcasecmp(params->ZBoundtype, "free") == 0) {
    params->ZBound = &Free_Boundary;
  } else if (strcasecmp(params->ZBoundtype, "periodic") == 0) {
    params->ZBound = &Periodic_Boundary;
  } else {
    fprintf(stderr, "Z boundary condition: %s not implemented \n", params->ZBoundtype);
    exit(2);
  }

  printf("xbound  %s\n", params->XBoundtype);
  printf("ybound  %s\n", params->YBoundtype);
  printf("zbound  %s\n", params->ZBoundtype);
  printf("\n");
}

int *Bulk_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  // Em geometria bulk, todos os pontos são do tipo 1 (interior)
  for (int ii = 0; ii < Nx * Ny * Nz; ii++) {
    pt[ii] = 1;
  }
  return pt;
}

float Bulk_Geometry::lattice_Potential(const nni fullni[7]) {
  float E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
  
  // Energia dos vizinhos primários
  E = Geometry::newman_neighbours(fullni);
  
  // Energia dos vizinhos secundários (se aplicável)
  if (params->neighbourKind > 1) {
    E += Geometry::second_neighbours(fullni);  // CORREÇÃO: second_neighbours
  }
  
  // Energia dos vizinhos terciários (se aplicável)
  if (params->neighbourKind == 3) {
    E += Geometry::third_neighbours(fullni);  // CORREÇÃO: third_neighbours
  }
  
  // Potencial elétrico (se aplicável)
  if (params->elecA != 0) {
    E += Electric_Potential(ni, params);
  }

  return E;
}

float Bulk_Geometry::Electric_Potential(float ni[3], Parameters& params) {
  // Implementação do potencial elétrico para geometria bulk
  return Potential::Electric_Potential(ni, &params);
}