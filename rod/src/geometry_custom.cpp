#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "../include/anchoring.h"
#include "../include/geometry.h"
#include "../include/parameters.h"
#include "../include/potential.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))

Custom_Geometry::Custom_Geometry(int *pt, Parameters *params) : Geometry(params) {
  printf("Geometry: Custom\n");

  ns = std::make_unique<float[]>(Nx * Ny * Nz * 3);  // CORREÇÃO: smart pointer
  pt = set_point_type_normals(pt, params);
  surfaces = std::vector<std::unique_ptr<Anchoring>>(nSurfaces);  // CORREÇÃO: smart pointers

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

int *Custom_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  int i, nn;
  int max_point_kind = 0, line_number, read_tester = 1;
  int NoV, pos_nx = -1, pos_ny = -1, pos_nz = -1, pos_pt = -1;
  
  // Abre arquivo de definição de geometria
  FILE *bound_input = fopen(params->bound_file_name, "r");
  if (bound_input == 0) {
    perror(params->bound_file_name);
    exit(2);
  }
  
  printf("Reading boundaries from %s\n", params->bound_file_name);
  
  // Lê cabeçalho para determinar número de colunas
  char line[500];  // CORREÇÃO: array na stack em vez de malloc
  std::string testline;
  if (fgets(line, 500, bound_input) == NULL) {
    fprintf(stderr, "Error reading header from boundary file\n");
    exit(2);
  }
  
  testline = line;
  NoV = 1 + ((int)std::count(testline.begin(), testline.end(), ','));
  std::unique_ptr<float[]> var = std::make_unique<float[]>(NoV);  // CORREÇÃO: smart pointer
  char strvar[500];
  
  // Volta ao início e identifica colunas
  rewind(bound_input);
  for (i = 0; i < NoV; i++) {
    char names[500];
    int dummy = fscanf(bound_input, "%[^,\n ],", names);
    if (strcasecmp(names, "nx") == 0)
      pos_nx = i;
    if (strcasecmp(names, "ny") == 0)
      pos_ny = i;
    if (strcasecmp(names, "nz") == 0)
      pos_nz = i;
    if (strcasecmp(names, "pt") == 0)
      pos_pt = i;
  }
  
  // Verifica se todas as colunas necessárias foram encontradas
  if (pos_nx == -1 || pos_ny == -1 || pos_nz == -1 || pos_pt == -1) {
    fprintf(stderr, "Missing required columns in boundary file\n");
    exit(2);
  }
  
  printf("Using the positions of nx %d ny %d nz %d and pt %d.\n", pos_nx, pos_ny, pos_nz, pos_pt);
  
  // Pula linha do cabeçalho
  if (fgets(line, 500, bound_input) == NULL) {
    fprintf(stderr, "Error reading after header from boundary file\n");
    exit(2);
  }
  
  // Inicializa todos os pontos como tipo 1 (bulk)
  for (i = 0; i < Nx * Ny * Nz; i++) {
    pt[i] = 1;
  }
  
  // Lê dados dos pontos de superfície
  line_number = 0;
  while (read_tester != EOF) {
    for (nn = 0; nn < NoV; nn++) {
      read_tester = fscanf(bound_input, "%[^,\n ],", strvar);
      if (read_tester == 1) {
        var[nn] = atof(strvar);
      } else {
        continue;
      }
    }
    
    // Verifica se leu dados válidos
    if (read_tester == EOF) break;
    
    int ii = (int)var[0];
    int jj = (int)var[1];
    int kk = (int)var[2];
    line_number++;
    
    // Verifica limites
    if (ii >= Nx || jj >= Ny || kk >= Nz) {
      fprintf(stderr, "Point out of the simulation box at line %d\n", line_number);
      exit(10);
    }
    
    // Verifica tipo de ponto válido
    if ((int)var[pos_pt] < 0) {
      fprintf(stderr,
              "Negative Point Type found at line %d\n\
Please, use 0 as empty space, 1 as bulk and 2+ as surface point type.\n",
              line_number);
      exit(10);
    }
    
    // Lê resto da linha
    if (fgets(line, 500, bound_input) == NULL && !feof(bound_input)) {
      fprintf(stderr, "Error reading line %d\n", line_number);
      exit(2);
    }
    
    // Atribui normais e tipo de ponto
    int index = ((Ny * kk + jj) * Nx + ii);
    ns[index * 3 + 0] = var[pos_nx];
    ns[index * 3 + 1] = var[pos_ny];
    ns[index * 3 + 2] = var[pos_nz];
    pt[index] = (int)var[pos_pt];
    max_point_kind = MAX(max_point_kind, pt[index]);
  }
  
  nSurfaces = max_point_kind - 1;
  printf("%d boundary(boundaries) found in %s!!\n\n", nSurfaces, params->bound_file_name);
  fflush(stdout);
  fclose(bound_input);
  
  return pt;
}

float Custom_Geometry::lattice_Potential(const nni fullni[7]) {
  float E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
  
  // Energia dos vizinhos primários
  E = Geometry::newman_neighbours(fullni);
  
  // Energia dos vizinhos terciários (se aplicável)
  if (params->neighbourKind == 3) {
    E += Geometry::third_neighbours(fullni);  // CORREÇÃO: third_neighbours
  }
  
  // Energia dos vizinhos secundários (se aplicável)
  if (params->neighbourKind > 1) {
    E += Geometry::second_neighbours(fullni);  // CORREÇÃO: second_neighbours
  }
  
  // Energia de superfície (se ponto de superfície)
  if (fullni[0].pt > 1) {
    float s[3] = {fullni[7].x, fullni[7].y, fullni[7].z};
    E += surfaces[fullni[0].pt - 2]->surface_potential(ni, s);
  }
  
  // Potencial elétrico (se aplicável)
  if (params->elecA != 0) {
    E += Electric_Potential(ni, params);
  }

  return E;
}

float Custom_Geometry::Electric_Potential(float ni[3], Parameters& params) {
  // Implementação do potencial elétrico para geometria customizada
  return Potential::Electric_Potential(ni, &params);
}