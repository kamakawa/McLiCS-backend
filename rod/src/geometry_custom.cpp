#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <fstream>  // ✅ Para std::ifstream

#include "../include/anchoring.h"
#include "../include/geometry.h"
#include "../include/parameters.h"
#include "../include/potential.h"

inline int MAX(int a, int b) { return (a > b) ? a : b; }

Custom_Geometry::Custom_Geometry(int *pt, Parameters *params) : Geometry(params) {
  printf("Geometry: Custom\n");

  ns = std::make_unique<float[]>(Nx * Ny * Nz * 3);
  for (int i = 0; i < Nx * Ny * Nz * 3; i++) {
    ns[i] = 0.0f;
  }

  pt = set_point_type_normals(pt, params);
  
  surfaces.clear();

  if (params->XBoundtype == "free")
    params->XBound = &Free_Boundary;
  else if (params->XBoundtype == "periodic")
    params->XBound = &Periodic_Boundary;
  else {
    fprintf(stderr, "X boundary condition: %s not implemented \n", params->XBoundtype.c_str());
    exit(2);
  }

  if (params->YBoundtype == "free")
    params->YBound = &Free_Boundary;
  else if (params->YBoundtype == "periodic")
    params->YBound = &Periodic_Boundary;
  else {
    fprintf(stderr, "Y boundary condition: %s not implemented \n", params->YBoundtype.c_str());
    exit(2);
  }

  if (params->ZBoundtype == "free")
    params->ZBound = &Free_Boundary;
  else if (params->ZBoundtype == "periodic")
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
  std::ifstream bound_input(params->bound_file_name);
  if (!bound_input.is_open()) {
    fprintf(stderr, "Cannot open boundary file: %s\n", params->bound_file_name.c_str());
    exit(2);
  }
  
  printf("Reading boundaries from %s\n", params->bound_file_name.c_str());
  
  std::string line;
  int NoV = 0;
  int pos_nx = -1, pos_ny = -1, pos_nz = -1, pos_pt = -1;
  
  if (std::getline(bound_input, line)) {
    NoV = 1 + static_cast<int>(std::count(line.begin(), line.end(), ','));
    
    size_t start = 0, end = 0;
    int col_index = 0;
    
    while (end != std::string::npos) {
      end = line.find(',', start);
      std::string column_name = line.substr(start, end - start);
      
      if (column_name == "nx") pos_nx = col_index;
      else if (column_name == "ny") pos_ny = col_index;
      else if (column_name == "nz") pos_nz = col_index;
      else if (column_name == "pt") pos_pt = col_index;
      
      start = (end == std::string::npos) ? end : end + 1;
      col_index++;
    }
  }
  
  if (pos_nx == -1 || pos_ny == -1 || pos_nz == -1 || pos_pt == -1) {
    fprintf(stderr, "Missing required columns in boundary file\n");
    exit(10);
  }
  
  printf("Using the positions of nx %d ny %d nz %d and pt %d.\n", pos_nx, pos_ny, pos_nz, pos_pt);
  
  for (int i = 0; i < Nx * Ny * Nz; i++) {
    pt[i] = 0; 
  }
  
  int line_number = 1; 
  int max_point_kind = 0;
  
  while (std::getline(bound_input, line)) {
    if (line.empty()) continue;
    
    std::vector<float> var;
    size_t start = 0, end = 0;
    
    while (end != std::string::npos) {
      end = line.find(',', start);
      std::string value_str = line.substr(start, end - start);
      if (!value_str.empty()) {
        var.push_back(std::stof(value_str));
      }
      start = (end == std::string::npos) ? end : end + 1;
    }
    
    if (var.size() != NoV) {
      fprintf(stderr, "Invalid number of columns at line %d\n", line_number);
      continue;
    }
    
    int ii = static_cast<int>(var[0]);
    int jj = static_cast<int>(var[1]);
    int kk = static_cast<int>(var[2]);
    
    if (ii >= Nx || jj >= Ny || kk >= Nz) {
      fprintf(stderr, "Point out of the simulation box at line %d\n", line_number);
      exit(10);
    }
    
    int point_type = static_cast<int>(var[pos_pt]);
    if (point_type < 0) {
      fprintf(stderr, 
              "Negative Point Type found at line %d\n"
              "Please, use 0 as empty space, 1 as bulk and 2+ as surface point type.\n",
              line_number);
      exit(10);
    }
    
    int idx = ((Ny * kk + jj) * Nx + ii);
    ns[idx * 3 + 0] = var[pos_nx];
    ns[idx * 3 + 1] = var[pos_ny];
    ns[idx * 3 + 2] = var[pos_nz];
    pt[idx] = point_type;
    
    max_point_kind = MAX(max_point_kind, point_type);
    line_number++;
  }
  
  nSurfaces = max_point_kind - 1;
  printf("%d boundary(boundaries) found in %s!!\n\n", nSurfaces, params->bound_file_name.c_str());
  fflush(stdout);
  
  Boundary_Init(params);
  
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
  
  if (fullni[0].pt > 1) {
    float s[3] = {ns[0], ns[1], ns[2]}; 
    E += surfaces[fullni[0].pt - 2]->surface_potential(ni, s);
  }
  
  if (params->elecA != 0) 
    E += Potential::Electric_Potential(ni, *params);

  return E;
}