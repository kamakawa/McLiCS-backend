#include "../include/geometry.h"

#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

#include "../include/anchoring.h"
#include "../include/parameters.h"
#include "../include/potential.h"

Geometry::Geometry(Parameters *params) : Nx(params->Nx), Ny(params->Ny), Nz(params->Nz) {
  this->params = params;
}

float Geometry::newman_neighbours(const nni fullni[]) {
  float rij[3];
  double E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};

  if (fullni[1].pt) {
    float nj[3] = {fullni[1].x, fullni[1].y, fullni[1].z};
    rij[0] = 1;
    rij[1] = 0;
    rij[2] = 0;
    E += bulk_potential(ni, nj, params, rij, 1);
  }

  if (fullni[2].pt) {
    float nj[3] = {fullni[2].x, fullni[2].y, fullni[2].z};
    rij[0] = -1;
    rij[1] = 0;
    rij[2] = 0;
    E += bulk_potential(ni, nj, params, rij, 1);
  }

  if (fullni[3].pt) {
    float nj[3] = {fullni[3].x, fullni[3].y, fullni[3].z};
    rij[0] = 0;
    rij[1] = 1;
    rij[2] = 0;
    E += bulk_potential(ni, nj, params, rij, 1);
  }

  if (fullni[4].pt) {
    float nj[3] = {fullni[4].x, fullni[4].y, fullni[4].z};
    rij[0] = 0;
    rij[1] = -1;
    rij[2] = 0;
    E += bulk_potential(ni, nj, params, rij, 1);
  }

  if (fullni[5].pt) {
    float nj[3] = {fullni[5].x, fullni[5].y, fullni[5].z};
    rij[0] = 0;
    rij[1] = 0;
    rij[2] = 1;
    E += bulk_potential(ni, nj, params, rij, 1);
  }

  if (fullni[6].pt) {
    float nj[3] = {fullni[6].x, fullni[6].y, fullni[6].z};
    rij[0] = 0;
    rij[1] = 0;
    rij[2] = -1;
    E += bulk_potential(ni, nj, params, rij, 1);
  }
  return E;
}

float Geometry::second_nerghbours(const nni fullni[]) {
  float rij[3];
  double E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
  const float isqrt2 = 0.707106781;
  
  for (int neighbor = 8; neighbor <= 19; neighbor++) {
    if (fullni[neighbor].pt) {
      float nj[3] = {fullni[neighbor].x, fullni[neighbor].y, fullni[neighbor].z};
      
      switch(neighbor) {
        case 8:  rij[0] = isqrt2; rij[1] = isqrt2; rij[2] = 0; break;
        case 9:  rij[0] = isqrt2; rij[1] = -isqrt2; rij[2] = 0; break;
        case 10: rij[0] = isqrt2; rij[1] = 0; rij[2] = isqrt2; break;
        case 11: rij[0] = isqrt2; rij[1] = 0; rij[2] = -isqrt2; break;
        case 12: rij[0] = -isqrt2; rij[1] = isqrt2; rij[2] = 0; break;
        case 13: rij[0] = -isqrt2; rij[1] = -isqrt2; rij[2] = 0; break;
        case 14: rij[0] = -isqrt2; rij[1] = 0; rij[2] = isqrt2; break;
        case 15: rij[0] = -isqrt2; rij[1] = 0; rij[2] = -isqrt2; break;
        case 16: rij[0] = 0; rij[1] = isqrt2; rij[2] = isqrt2; break;
        case 17: rij[0] = 0; rij[1] = isqrt2; rij[2] = -isqrt2; break;
        case 18: rij[0] = 0; rij[1] = -isqrt2; rij[2] = isqrt2; break;
        case 19: rij[0] = 0; rij[1] = -isqrt2; rij[2] = -isqrt2; break;
      }
      
      E += bulk_potential(ni, nj, params, rij, 2);
    }
  }
  return E;
}

float Geometry::third_nerghbours(const nni fullni[]) {
  float rij[3];
  double E = 0;
  float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
  const float isqrt3 = 0.577350269;
  
  for (int neighbor = 20; neighbor <= 27; neighbor++) {
    if (fullni[neighbor].pt) {
      float nj[3] = {fullni[neighbor].x, fullni[neighbor].y, fullni[neighbor].z};
      
      switch(neighbor) {
        case 20: rij[0] = isqrt3; rij[1] = isqrt3; rij[2] = isqrt3; break;
        case 21: rij[0] = isqrt3; rij[1] = isqrt3; rij[2] = -isqrt3; break;
        case 22: rij[0] = isqrt3; rij[1] = -isqrt3; rij[2] = isqrt3; break;
        case 23: rij[0] = isqrt3; rij[1] = -isqrt3; rij[2] = -isqrt3; break;
        case 24: rij[0] = -isqrt3; rij[1] = isqrt3; rij[2] = isqrt3; break;
        case 25: rij[0] = -isqrt3; rij[1] = isqrt3; rij[2] = -isqrt3; break;
        case 26: rij[0] = -isqrt3; rij[1] = -isqrt3; rij[2] = isqrt3; break;
        case 27: rij[0] = -isqrt3; rij[1] = -isqrt3; rij[2] = -isqrt3; break;
      }
      
      E += bulk_potential(ni, nj, params, rij, 3);
    }
  }
  return E;
}

void Geometry::Boundary_Init(Parameters *params) {
  std::string anchoring;
  for (int ii = 0; ii < nSurfaces; ii++) {
    try {
      anchoring = params->anchoring_type.at(ii);
    } catch (std::out_of_range dummy_var) {
      std::cout << "You must define " << nSurfaces << " boundaries.\nPlease review your input file.\nAborting the program.\n\n";
      exit(0);
    }
    
    if (strcasecmp(anchoring.c_str(), "rp") == 0)
      surfaces.push_back(std::make_unique<RP_Anchoring>(params, ii));
    else if (strcasecmp(anchoring.c_str(), "fg") == 0)
      surfaces.push_back(std::make_unique<FG_Anchoring>(params, ii));
    else if (strcasecmp(anchoring.c_str(), "homeotropic") == 0)
      surfaces.push_back(std::make_unique<Homeotropic_Anchoring>(params, ii));
    else if (strcasecmp(anchoring.c_str(), "strong") == 0)
      surfaces.push_back(std::make_unique<Strong_Anchoring>(params, ii));
    else if (strcasecmp(anchoring.c_str(), "rp_ghrl") == 0)
      surfaces.push_back(std::make_unique<RP_Anchoring_GHRL>(params, ii));
    else if (strcasecmp(anchoring.c_str(), "fg_ghrl") == 0)
      surfaces.push_back(std::make_unique<FG_Anchoring_GHRL>(params, ii));
    else if (strcasecmp(anchoring.c_str(), "homeotropic_ghrl") == 0)
      surfaces.push_back(std::make_unique<Homeotropic_Anchoring_GHRL>(params, ii));
    else if (strcasecmp(anchoring.c_str(), "strong_ghrl") == 0)
      surfaces.push_back(std::make_unique<Strong_Anchoring_GHRL>(params, ii));
    else {
      printf("%s boundary condition is not defined\n", anchoring.c_str());
      exit(2);
    }
  }
}


Bulk_Geometry::Bulk_Geometry(int *pt, Parameters *params) : Geometry(params) {
  set_point_type_normals(pt, params);
}

int* Bulk_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  int total_size = params->Nx * params->Ny * params->Nz;
  for (int i = 0; i < total_size; i++) {
    pt[i] = 0; 
  }
  return pt;
}

float Bulk_Geometry::lattice_Potential(const nni ni[7]) {
  return newman_neighbours(ni);
}

Slab_Geometry::Slab_Geometry(int *pt, Parameters *params) : Geometry(params) {
  set_point_type_normals(pt, params);
}

int* Slab_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  int total_size = params->Nx * params->Ny * params->Nz;
  
  for (int k = 0; k < params->Nz; k++) {
    for (int j = 0; j < params->Ny; j++) {
      for (int i = 0; i < params->Nx; i++) {
        int idx = i + j * params->Nx + k * params->Nx * params->Ny;
        
        if (k == 0 || k == params->Nz - 1) {
          pt[idx] = 1; 
        } else {
          pt[idx] = 0; 
        }
      }
    }
  }
  return pt;
}

float Slab_Geometry::lattice_Potential(const nni ni[7]) {
  float energy = newman_neighbours(ni);
  
  if (ni[0].pt == 1) {
    for (auto& surface : surfaces) {
      float s[3] = {0, 0, 1}; 
      float ni_vec[3] = {ni[0].x, ni[0].y, ni[0].z};
      energy += surface->surface_potential(ni_vec, s);
    }
  }
  
  return energy;
}

Sphere_Geometry::Sphere_Geometry(int *pt, Parameters *params) : Geometry(params) {
  set_point_type_normals(pt, params);
}

int* Sphere_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  int total_size = params->Nx * params->Ny * params->Nz;
  float center_x = params->Nx / 2.0f;
  float center_y = params->Ny / 2.0f; 
  float center_z = params->Nz / 2.0f;
  float radius = std::min({center_x, center_y, center_z}) - 1;
  
  for (int k = 0; k < params->Nz; k++) {
    for (int j = 0; j < params->Ny; j++) {
      for (int i = 0; i < params->Nx; i++) {
        int idx = i + j * params->Nx + k * params->Nx * params->Ny;
        float dx = i - center_x;
        float dy = j - center_y;
        float dz = k - center_z;
        float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        if (std::abs(dist - radius) < 1.0f) {
          pt[idx] = 1; 
        } else if (dist < radius) {
          pt[idx] = 0; 
        } else {
          pt[idx] = -1; 
        }
      }
    }
  }
  return pt;
}

float Sphere_Geometry::lattice_Potential(const nni ni[7]) {
  float energy = newman_neighbours(ni);
  
  if (ni[0].pt == 1) { 
    float s[3] = {0, 0, 1}; 
    float ni_vec[3] = {ni[0].x, ni[0].y, ni[0].z};
    
    for (auto& surface : surfaces) {
      energy += surface->surface_potential(ni_vec, s);
    }
  }
  
  return energy;
}

Custom_Geometry::Custom_Geometry(int *pt, Parameters *params) : Geometry(params) {
  set_point_type_normals(pt, params);
}

int* Custom_Geometry::set_point_type_normals(int *pt, Parameters *params) {
  int total_size = params->Nx * params->Ny * params->Nz;
  
  for (int i = 0; i < total_size; i++) {
    pt[i] = 0; 
  }
  
  return pt;
}

float Custom_Geometry::lattice_Potential(const nni ni[7]) {
  float energy = 0.0f;
  
  switch(params->neighbourKind) {
    case 1:
      energy = newman_neighbours(ni);
      break;
    case 2:
      energy = newman_neighbours(ni) + second_nerghbours(ni);
      break;
    case 3:
      energy = newman_neighbours(ni) + second_nerghbours(ni) + third_nerghbours(ni);
      break;
    default:
      energy = newman_neighbours(ni);
  }
  
  return energy;
}