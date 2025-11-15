#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <memory> 

#include "../include/anchoring.h"
#include "../include/geometry.h"
#include "../include/parameters.h"
#include "../include/potential.h"

#define MAX(a, b) (a > b ? a : b)

Custom_Geometry::Custom_Geometry(int *pt, Parameters *params) : Geometry(params) {
  printf("Geometry: Custom\n");

  ns = std::unique_ptr<float[]>(new float[Nx * Ny * Nz * 3]()); 
  
  pt = set_point_type_normals(pt, params); 
  
  if (strcasecmp(params->XBoundtype.c_str(), "free") == 0)
    params->XBound = &Free_Boundary;
  else if (strcasecmp(params->XBoundtype.c_str(), "periodic") == 0)
    params->XBound = &Periodic_Boundary;
  else {
    fprintf(stderr, "X boundary condition: %s not implemented \n", params->XBoundtype.c_str());
    exit(2);
  }

  if (strcasecmp(params->YBoundtype.c_str(), "free") == 0)
    params->YBound = &Free_Boundary;
  else if (strcasecmp(params->YBoundtype.c_str(), "periodic") == 0)
    params->YBound = &Periodic_Boundary;
  else {
    fprintf(stderr, "Y boundary condition: %s not implemented \n", params->YBoundtype.c_str());
    exit(2);
  }

  if (strcasecmp(params->ZBoundtype.c_str(), "free") == 0)
    params->ZBound = &Free_Boundary;
  else if (strcasecmp(params->ZBoundtype.c_str(), "periodic") == 0)
    params->ZBound = &Periodic_Boundary;
  else {
    fprintf(stderr, "Z boundary condition: %s not implemented \n", params->ZBoundtype.c_str());
    exit(2);
  }

  printf("xbound  %s\n", params->XBoundtype.c_str());
  printf("ybound  %s\n", params->YBoundtype.c_str());
  printf("zbound  %s\n", params->ZBoundtype.c_str());
  printf("\n");
}

int *Custom_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  int i, nn;
  int max_point_kind = 0, line_number, read_tester = 1;
  int NoV, pos_nx, pos_ny, pos_nz, pos_pt;
  
  FILE *bound_input = fopen(params->bound_file_name.c_str(), "r");
  if (bound_input == 0) {
    perror(params->bound_file_name.c_str());
    exit(2);
  }
  printf("Reading boundaries from %s\n", params->bound_file_name.c_str());

  std::unique_ptr<char[]> line_ptr(new char[500]); 
  char* line = line_ptr.get();
  
  std::string testline;
  line = fgets(line, 500, bound_input);
  testline = line;
  NoV = 1 + ((int) std::count(testline.begin(), testline.end(), ','));
  
  std::unique_ptr<float[]> var_ptr(new float[NoV]);
  float* var = var_ptr.get();
  
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
  for (i = 0; i < Nx * Ny * Nz; i++) {
    pt[i] = 1;
  }
  
  float *ns_ptr = ns.get(); 

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
    if (ii >= Nx || jj >= Ny || kk >= Nz) {
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
    ns_ptr[((Ny * kk + jj) * Nx + ii) * 3 + 0] = var[pos_nx];
    ns_ptr[((Ny * kk + jj) * Nx + ii) * 3 + 1] = var[pos_ny];
    ns_ptr[((Ny * kk + jj) * Nx + ii) * 3 + 2] = var[pos_nz];
    pt[(Ny * kk + jj) * Nx + ii] = (int)var[pos_pt];
    max_point_kind = MAX(max_point_kind, pt[(Ny * kk + jj) * Nx + ii]);
  }
  nSurfaces = max_point_kind - 1;
  printf("%d boundary(boundaries) found in %s!!\n\n", nSurfaces, params->bound_file_name.c_str());
  fflush(stdout);
  fclose(bound_input);
  
  return pt;
}

float Custom_Geometry::lattice_Potential(const nni fullni[7]) {
  float E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
  
  E = Geometry::newman_neighbours(fullni);

  if (params->neighbourKind == 3)
    E += Geometry::third_nerghbours(fullni);
  if (params->neighbourKind > 1)
    E += Geometry::second_nerghbours(fullni);
    
  if (fullni[0].pt > 1){
    float s[3] = {fullni[7].x, fullni[7].y, fullni[7].z};
    
    E += surfaces[fullni[0].pt - 2]->surface_potential(ni, s);
  }
  
  if (params->elecA!=0) 
    E += Potential::Electric_Potential(ni, *params);

  return E;
}