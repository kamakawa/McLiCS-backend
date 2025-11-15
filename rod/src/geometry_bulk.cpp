#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <memory>

#include "../include/anchoring.h"
#include "../include/geometry.h"
#include "../include/parameters.h"
#include "../include/potential.h"

Bulk_Geometry::Bulk_Geometry(int *pt, Parameters *params) : Geometry(params) {
  printf("Geometry: Bulk\n");
  
  ns = std::unique_ptr<float[]>(new float[Nx * Ny * Nz * 3]()); 

  nSurfaces = 0;

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

int *Bulk_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  for (int ii = 0; ii < Nx * Ny * Nz; ii++) pt[ii] = 1;

  return pt;
}

float Bulk_Geometry::lattice_Potential(const nni fullni[7]) { 
  float E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
  
  E = Geometry::newman_neighbours(fullni);
  
  if (params->neighbourKind > 1)
    E += Geometry::second_nerghbours(fullni);
  if (params->neighbourKind == 3)
    E += Geometry::third_nerghbours(fullni);
    
  if (params->elecA!=0) 
    E += Potential::Electric_Potential(ni, *params);

  return E;
}
