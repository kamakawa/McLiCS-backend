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
  ns = std::make_unique<float[]>(Nx * Ny * Nz * 3);  // CORREÇÃO: smart pointer
  surfaces = std::vector<std::unique_ptr<Anchoring>>(nSurfaces);  // CORREÇÃO: smart pointers
  pt = set_point_type_normals(pt, params);

  // Força condições de contorno free na direção Z para geometria slab
  strcpy(params->ZBoundtype, "free");  // CORREÇÃO: strcpy em vez de sprintf
  
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

  printf("xbound  %s\n", params->XBoundtype);
  printf("ybound  %s\n", params->YBoundtype);
  printf("zbound  %s (forced free for slab geometry)\n", params->ZBoundtype);
  printf("\n");
}

int *Slab_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  // Configura pontos da superfície inferior (z = 0)
  for (int ii = 0; ii < Nx; ii++) {
    for (int jj = 0; jj < Ny; jj++) {
      int kk = 0;
      // Superfície inferior: tipo 2, normal para baixo (0,0,-1)
      int index = ii + Nx * (jj + Ny * kk);
      pt[index] = 2;
      ns[index * 3 + 0] = 0;
      ns[index * 3 + 1] = 0;
      ns[index * 3 + 2] = -1;
    }
  }
  
  // Configura pontos do bulk (1 < z < Nz-1)
  for (int ii = 0; ii < Nx; ii++) {
    for (int jj = 0; jj < Ny; jj++) {
      for (int kk = 1; kk < Nz - 1; kk++) {
        int index = ii + Nx * (jj + Ny * kk);
        pt[index] = 1;  // Bulk
      }
    }
  }
  
  // Configura pontos da superfície superior (z = Nz-1)
  for (int ii = 0; ii < Nx; ii++) {
    for (int jj = 0; jj < Ny; jj++) {
      int kk = Nz - 1;
      // Superfície superior: tipo 3, normal para cima (0,0,1)
      int index = ii + Nx * (jj + Ny * kk);
      pt[index] = 3;
      ns[index * 3 + 0] = 0;
      ns[index * 3 + 1] = 0;
      ns[index * 3 + 2] = 1;
    }
  }
  
  return pt;
}

float Slab_Geometry::lattice_Potential(const nni fullni[7]) {
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
  
  // Energia de superfície (se ponto de superfície)
  float s[3] = {fullni[7].x, fullni[7].y, fullni[7].z};
  if (fullni[0].pt > 1) {
    E += surfaces[fullni[0].pt - 2]->surface_potential(ni, s);
  }
  
  // Potencial elétrico (se aplicável)
  if (params->elecA != 0) {
    E += Electric_Potential(ni, params);
  }

  return E;
}

float Slab_Geometry::Electric_Potential(float ni[3], Parameters& params) {
  // Implementação do potencial elétrico para geometria slab
  return Potential::Electric_Potential(ni, &params);
}