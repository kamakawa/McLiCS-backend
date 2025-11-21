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

// Construtor da Geometria Esferica
Sphere_Geometry::Sphere_Geometry(int *pt, Parameters *params) : Geometry(params) {
  printf("Geometry: Sphere\n");
  
  nSurfaces = 1; // Apenas uma superficie (a casca da esfera)
  
  // Forca condicoes de contorno livres em todas as direcoes
  sprintf(params->XBoundtype, "free");
  sprintf(params->YBoundtype, "free");
  sprintf(params->ZBoundtype, "free");
  
  ns = (float *)calloc(Nx * Ny * Nz * 3, sizeof(float));
  surfaces = std::vector<class Anchoring *>(nSurfaces);
  
  pt = set_point_type_normals(pt, params);
  printf("\n");
}

// Define a geometria esferica (Bulk, Borda e Vazio)
int *Sphere_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  float Rx, Ry, Rz;
  float Hx = Nx / 2;
  float Hy = Ny / 2;
  float Hz = Nz / 2; // Define o centro da caixa como centro da esfera

  for (int ii = 0; ii < Nx; ii++) {
    for (int jj = 0; jj < Ny; jj++) {
      for (int kk = 0; kk < Nz; kk++) {
        
        // Distancia do ponto atual ao centro
        Rx = ii - Hx;
        Ry = jj - Hy;
        Rz = kk - Hz;
        float radius = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);
        
        int idx = ii + Nx * (jj + Ny * kk);

        // --- Logica de Definicao da Esfera ---
        if (radius < Hz - 1) {
          // Parte interna (Bulk)
          pt[idx] = 1;
        } 
        else if (radius < Hz + 1) {
          // Casca (Superficie)
          pt[idx] = 2;
          ns[idx * 3 + 0] = Rx / radius;
          ns[idx * 3 + 1] = Ry / radius;
          ns[idx * 3 + 2] = Rz / radius;
        } 
        else {
          // Fora da esfera (Vazio)
          pt[idx] = 0;
        }
      }
    }
  }
  return pt;
}

// Calculo do Potencial
float Sphere_Geometry::latice_Potential(const nni fullni[7]) {
  float rij[3];
  float E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
  float s[3] = {fullni[7].x, fullni[7].y, fullni[7].z};

  E = Geometry::newman_neighbours(fullni);

  if (params->neighbourKind > 1)
    E += Geometry::second_nerghbours(fullni);
  
  if (params->neighbourKind == 3)
    E += Geometry::third_nerghbours(fullni);
  
  // Interacao com a superficie da esfera
  if (fullni[0].pt > 1)
    E += surfaces[fullni[0].pt - 2]->surface_potential(ni, s);
  
  if (params->elecA != 0) 
    E += Potential::Electric_Potential(ni, params);

  return E;
}