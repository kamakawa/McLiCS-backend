#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// --- Project Includes ---
#include "../include/anchoring.h"
#include "../include/geometry.h"
#include "../include/parameters.h"
#include "../include/potential.h"

// Construtor da Geometria Slab (Placa)
Slab_Geometry::Slab_Geometry(int *pt, Parameters *params) : Geometry(params) {
  printf("Geometry: Slab\n");
  
  nSurfaces = 2; 
  ns = (float *)calloc(Nx * Ny * Nz * 3, sizeof(float));
  surfaces = std::vector<class Anchoring *>(nSurfaces);
  
  pt = set_point_type_normals(pt, params);

  // Forca condicao de contorno Livre em Z (caracteristica do Slab)
  sprintf(params->ZBoundtype, "free");
  
  // --- Configuracao da Condicao de Contorno em X ---
  if (strcasecmp(params->XBoundtype, "free") == 0) {
    params->XBound = &Potential::Free_Boundary;
  } else if (strcasecmp(params->XBoundtype, "periodic") == 0) {
    params->XBound = &Potential::Periodic_Boundary;
  } else {
    fprintf(stderr, "X boundary condition: %s not implemented \n", params->XBoundtype);
    exit(2);
  }

  // --- Configuracao da Condicao de Contorno em Y ---
  if (strcasecmp(params->YBoundtype, "free") == 0) {
    params->YBound = &Potential::Free_Boundary;
  } else if (strcasecmp(params->YBoundtype, "periodic") == 0) {
    params->YBound = &Potential::Periodic_Boundary;
  } else {
    fprintf(stderr, "Y boundary condition: %s not implemented \n", params->YBoundtype);
    exit(2);
  }

  printf("xbound  %s\n", params->XBoundtype);
  printf("ybound  %s\n", params->YBoundtype);
  printf("\n");
}

// Configura os tipos de pontos e vetores normais para a placa
int *Slab_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  int kk;
  
  for (int ii = 0; ii < Nx; ii++) {
    for (int jj = 0; jj < Ny; jj++) {
      
      kk = 0;
      int idx_bottom = ii + Nx * (jj + Ny * kk);
      
      pt[idx_bottom] = 2; 
      ns[idx_bottom * 3 + 0] = 0;
      ns[idx_bottom * 3 + 1] = 0;
      ns[idx_bottom * 3 + 2] = -1;

      for (kk = 1; kk < Nz - 1; kk++) {
        pt[ii + Nx * (jj + Ny * kk)] = 1; 
      }

      kk = Nz - 1;
      int idx_top = ii + Nx * (jj + Ny * kk);
      
      pt[idx_top] = 3; 
      ns[idx_top * 3 + 0] = 0;
      ns[idx_top * 3 + 1] = 0;
      ns[idx_top * 3 + 2] = 1;
    }
  }
  return pt;
}

// Calculo do Potencial
float Slab_Geometry::latice_Potential(const nni fullni[7]) {
  float rij[3];
  float E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};

  E = Geometry::newman_neighbours(fullni);

  if (params->neighbourKind > 1)
    E += Geometry::second_nerghbours(fullni);
  
  if (params->neighbourKind == 3)
    E += Geometry::third_nerghbours(fullni);
  
  // Interacao de Superficie
  if (fullni[0].pt > 1) {
    float s[3] = {fullni[7].x, fullni[7].y, fullni[7].z};
    E += surfaces[fullni[0].pt - 2]->surface_potential(ni, s);
  }

  if (params->elecA != 0) 
    E += Potential::Electric_Potential(ni, params);

  return E;
}