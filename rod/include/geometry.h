#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <gsl/gsl_rng.h>
#include <iostream>
#include <memory>
#include <vector>

#include "../include/define.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

int Free_Boundary(int &ii, int NN);
int Periodic_Boundary(int &ii, int NN);

class Geometry {
 public:
   Geometry(Parameters *params) 
      : Nx(params->Nx), Ny(params->Ny), Nz(params->Nz),  
        nSurfaces(0), bulk_potential(nullptr), params(params) {}

  int Nx, Ny, Nz;
  std::vector<std::unique_ptr<class Anchoring>> surfaces;
  int nSurfaces;
  float newman_neighbours(const nni fullni[]);
  float second_nerghbours(const nni fullni[]);
  float third_nerghbours(const nni fullni[]);
  void Boundary_Init(Parameters *params);
  virtual int *set_point_type_normals(int *pt, Parameters *params) = 0;
  virtual float lattice_Potential(const nni ni[7]) = 0;
  
  float (*bulk_potential)(float ni[3], float nj[3], Parameters *params, float rij[3], int nk);
  std::unique_ptr<float[]> ns;
  Parameters *params;

  virtual ~Geometry() = default;
};

class Bulk_Geometry : public Geometry {
 public:
  Bulk_Geometry(int *pt, Parameters *params);
  float lattice_Potential(const nni ni[7]) override;

 private:
  int *set_point_type_normals(int *pt, Parameters *params) override;
};

class Slab_Geometry : public Geometry {
 public:
  Slab_Geometry(int *pt, Parameters *params);
  float lattice_Potential(const nni ni[7]) override;

 private:
  int *set_point_type_normals(int *pt, Parameters *params) override;
};

class Sphere_Geometry : public Geometry {
 public:
  Sphere_Geometry(int *pt, Parameters *params);
  float lattice_Potential(const nni ni[7]) override;

 private:
  int *set_point_type_normals(int *pt, Parameters *params) override;
};

class Custom_Geometry : public Geometry {
 public:
  Custom_Geometry(int *pt, Parameters *params);
  float lattice_Potential(const nni ni[7]) override;

 private:
  int *set_point_type_normals(int *pt, Parameters *params) override;
};

#endif