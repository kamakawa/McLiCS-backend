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
#define MAX(a, b) (a > b ? a : b)
Custom_Geometry::Custom_Geometry(int *pt, Parameters *params) : Geometry(params) {
  printf("Geometry: Custom\n");

  // Alteração 1, 2, 3: Acesso às dimensões da grade aninhadas em 'lattice'
  ns = (float *)calloc(params->lattice.Nx * params->lattice.Ny * params->lattice.Nz * 3, sizeof(float));
  pt = set_point_type_normals(pt, params);
  surfaces = std::vector<class Anchoring *>(nSurfaces);

  // Alteração 4, 5: Acesso a XBoundtype e XBound aninhados em 'lattice'
  if (strcasecmp(params->lattice.XBoundtype, "free") == 0)
    params->lattice.XBound = &Free_Boundary; // Alteração 6: Chamada de função do namespace Potential
  else if (strcasecmp(params->lattice.XBoundtype, "periodic") == 0)
    params->lattice.XBound = &Periodic_Boundary; // Alteração 7: Chamada de função do namespace Potential
  else {
    fprintf(stderr, "X boundary condition: %s not implemented \n", params->lattice.XBoundtype);
    exit(2);
  }

  // Alteração 8, 9: Acesso a YBoundtype e YBound aninhados em 'lattice'
  if (strcasecmp(params->lattice.YBoundtype, "free") == 0)
    params->lattice.YBound = &Free_Boundary; // Alteração 10: Chamada de função do namespace Potential
  else if (strcasecmp(params->lattice.YBoundtype, "periodic") == 0)
    params->lattice.YBound = &Periodic_Boundary; // Alteração 11: Chamada de função do namespace Potential
  else {
    fprintf(stderr, "Y boundary condition: %s not implemented \n", params->lattice.YBoundtype);
    exit(2);
  }

  // Alteração 12, 13: Acesso a ZBoundtype e ZBound aninhados em 'lattice'
  if (strcasecmp(params->lattice.ZBoundtype, "free") == 0)
    params->lattice.ZBound = &Free_Boundary; // Alteração 14: Chamada de função do namespace Potential
  else if (strcasecmp(params->lattice.ZBoundtype, "periodic") == 0)
    params->lattice.ZBound = &Periodic_Boundary; // Alteração 15: Chamada de função do namespace Potential
  else {
    fprintf(stderr, "Z boundary condition: %s not implemented \n", params->lattice.ZBoundtype);
    exit(2);
  }

  // Alteração 16, 17, 18: Acesso aos Boundarytypes aninhados em 'lattice'
  printf("xbound  %s\n", params->lattice.XBoundtype);
  printf("ybound  %s\n", params->lattice.YBoundtype);
  printf("zbound  %s\n", params->lattice.ZBoundtype);
  printf("\n");
}

int *Custom_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  int i, nn;
  int max_point_kind = 0, line_number, read_tester = 1;
  // Alteração 19, 20, 21: Acesso às dimensões da grade aninhadas em 'lattice'
  int NoV, pos_nx, pos_ny, pos_nz, pos_pt;
  FILE *bound_input = fopen(params->surface.bound_file_name, "r"); // Alteração 22: Acesso a bound_file_name aninhado em 'surface'
  if (bound_input == 0) {
    perror(params->surface.bound_file_name); // Alteração 23: Acesso a bound_file_name aninhado em 'surface'
    exit(2);
  }
  printf("Reading boundaries from %s\n", params->surface.bound_file_name); // Alteração 24: Acesso a bound_file_name aninhado em 'surface'
  char *line = (char *)malloc(500);
  std::string testline;
  line = fgets(line, 500, bound_input);
  testline = line;
  NoV = 1 + ((int) std::count(testline.begin(), testline.end(), ','));
  float *var = new float[NoV];
  char strvar[500];
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
  printf("Using the positions of nx %d ny %d nz %d and pt %d.\n ", pos_nx, pos_ny, pos_nz, pos_pt);
  line_number = 0;
  line = fgets(line, 500, bound_input);
  // Alteração 25, 26, 27: Acesso às dimensões da grade aninhadas em 'lattice'
  for (i = 0; i < params->lattice.Nx * params->lattice.Ny * params->lattice.Nz; i++) {
    pt[i] = 1;
  }
  while (read_tester != EOF) {
    for (nn = 0; nn < NoV; nn++) {
      read_tester = fscanf(bound_input, "%[^,\n ],", strvar);
      if (read_tester == 1)
        var[nn] = atof(strvar);
      else
        continue;
    }
    int ii = (int)var[0];
    int jj = (int)var[1];
    int kk = (int)var[2];
    line_number++;
    // Alteração 28, 29, 30: Acesso às dimensões da grade aninhadas em 'lattice'
    if (ii >= params->lattice.Nx || jj >= params->lattice.Ny || kk >= params->lattice.Nz) {
      fprintf(stderr, "Point out of the simulation box at line %d\n", line_number);
      exit(10);
    }
    if ((int)var[pos_pt] < 0) {
      fprintf(stderr,
              "Negative Point Type founde at line %d\n\
Please, use 0 as empty space, 1 as bulk and 2+ as surface point type.\n",
              line_number);
      exit(10);
    }
    line = fgets(line, 500, bound_input);
    // Alteração 31, 32, 33: Acesso às dimensões da grade aninhadas em 'lattice'
    ns[((params->lattice.Ny * kk + jj) * params->lattice.Nx + ii) * 3 + 0] = var[pos_nx];
    ns[((params->lattice.Ny * kk + jj) * params->lattice.Nx + ii) * 3 + 1] = var[pos_ny];
    ns[((params->lattice.Ny * kk + jj) * params->lattice.Nx + ii) * 3 + 2] = var[pos_nz];
    pt[(params->lattice.Ny * kk + jj) * params->lattice.Nx + ii] = (int)var[pos_pt];
    max_point_kind = MAX(max_point_kind, pt[(params->lattice.Ny * kk + jj) * params->lattice.Nx + ii]);
  }
  nSurfaces = max_point_kind - 1;
  printf("%d boundary(boundaries) found in %s!!\n\n", nSurfaces, params->surface.bound_file_name); // Alteração 34: Acesso a bound_file_name aninhado em 'surface'
  fflush(stdout);
  fclose(bound_input);
  return pt;
}

float Custom_Geometry::latice_Potential(const nni fullni[7]) {
  float rij[3];
  float E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
  
  E = Geometry::newman_neighbours(fullni);

  // Alteração 35: Acesso a neighbourKind aninhado em 'neighbourhood'
  if (params->neighbourhood.neighbourKind == 3)
    E += Geometry::third_nerghbours(fullni);
    
  // Alteração 36: Acesso a neighbourKind aninhado em 'neighbourhood'
  if (params->neighbourhood.neighbourKind > 1)
    E += Geometry::second_nerghbours(fullni);
    
  if (fullni[0].pt > 1){
    float s[3] = {fullni[7].x, fullni[7].y, fullni[7].z};
    E += surfaces[fullni[0].pt - 2]->surface_potential(ni, s);
  }
  
  // Alteração 37: Acesso a elecA aninhado em 'electric'
  // Alteração 38: Chamada de função do namespace Potential
  if (params->electric.elecA!=0) 
    E+=Potential::Electric_Potential(ni,params);

  return E;
}