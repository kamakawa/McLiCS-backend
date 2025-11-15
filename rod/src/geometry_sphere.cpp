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
  
  // Força condições de contorno free em todas as direções para geometria esférica
  strcpy(params->XBoundtype, "free");  // CORREÇÃO: strcpy em vez de sprintf
  strcpy(params->YBoundtype, "free");
  strcpy(params->ZBoundtype, "free");
  
  ns = std::make_unique<float[]>(Nx * Ny * Nz * 3);  // CORREÇÃO: smart pointer
  surfaces = std::vector<std::unique_ptr<Anchoring>>(nSurfaces);  // CORREÇÃO: smart pointers
  pt = set_point_type_normals(pt, params);
  
  printf("xbound  %s (forced free for sphere geometry)\n", params->XBoundtype);
  printf("ybound  %s (forced free for sphere geometry)\n", params->YBoundtype);
  printf("zbound  %s (forced free for sphere geometry)\n", params->ZBoundtype);
  printf("\n");
}

int *Sphere_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  // Centro da esfera
  const float Hx = Nx / 2.0f;
  const float Hy = Ny / 2.0f; 
  const float Hz = Nz / 2.0f;
  
  // Raio da esfera (usa a menor dimensão como referência)
  const float sphere_radius = std::min({Hx, Hy, Hz}) - 1.0f;
  
  for (int ii = 0; ii < Nx; ii++) {
    for (int jj = 0; jj < Ny; jj++) {
      for (int kk = 0; kk < Nz; kk++) {
        // Vetor do centro para o ponto atual
        float Rx = ii - Hx;
        float Ry = jj - Hy;
        float Rz = kk - Hz;
        float radius = sqrtf(Rx * Rx + Ry * Ry + Rz * Rz);
        
        int index = ii + Nx * (jj + Ny * kk);
        
        if (radius < sphere_radius - 1.0f) {
          // Interior da esfera: bulk
          pt[index] = 1;
        } else if (radius < sphere_radius + 1.0f) {
          // Superfície da esfera: tipo 2 com normal radial
          pt[index] = 2;
          ns[index * 3 + 0] = Rx / radius;  // Normal x
          ns[index * 3 + 1] = Ry / radius;  // Normal y  
          ns[index * 3 + 2] = Rz / radius;  // Normal z
        } else {
          // Exterior da esfera: vazio
          pt[index] = 0;
        }
      }
    }
  }
  return pt;
}

float Sphere_Geometry::lattice_Potential(const nni fullni[7]) {
  float E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
  float s[3] = {fullni[7].x, fullni[7].y, fullni[7].z};
  
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
  if (fullni[0].pt > 1) {
    E += surfaces[fullni[0].pt - 2]->surface_potential(ni, s);
  }
  
  // Potencial elétrico (se aplicável)
  if (params->elecA != 0) {
    E += Electric_Potential(ni, params);
  }

  return E;
}

float Sphere_Geometry::Electric_Potential(float ni[3], Parameters& params) {
  // Implementação do potencial elétrico para geometria esférica
  return Potential::Electric_Potential(ni, &params);
}