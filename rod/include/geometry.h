//Modifiquei
//Adicionado destrutor virtual em Geometry
//Todos os membros inicializados explicitamente com valores seguros
//Elimina valores "lixo" e comportamentos indefinidos
#ifndef GEOMETRY_H_
#define GEOMETRY_H_
#include <gsl/gsl_rng.h>

#include <iostream>
#include <vector>

#include "../include/define.h"
//~ #include "../include/evolve.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

namespace Potential { 
  int Free_Boundary(int &ii, int NN); 
  int Periodic_Boundary(int &ii, int NN);
  float Electric_Potential(const float ni[3], const Parameters* params);
}
class Geometry {
 public:
   Geometry(Parameters *params) 
      : Nx(params->lattice.Nx), Ny(params->lattice.Ny), Nz(params->lattice.Nz), // Inicialização explícita
        nSurfaces(0), surfaces(), bulk_potential(nullptr),
        ns(nullptr), params(params) {
  }

  int Nx, Ny, Nz;
  std::vector<class Anchoring *> surfaces;
  int nSurfaces;
  float newman_neighbours(const nni fullni[]);
  float second_nerghbours(const nni fullni[]);
  float third_nerghbours(const nni fullni[]);
  void Boundary_Init(Parameters *params);
  virtual int *set_point_type_normals(int *pt, Parameters *params) = 0;
  virtual float latice_Potential(const nni ni[7]) = 0;
  float (*bulk_potential)(const float ni[3], const float nj[3], const Parameters *params, const float rij[3], int nk);  float *ns;
  Parameters *params;

  ~Geometry() = default;
};

class Bulk_Geometry : public Geometry {
 public:
  Bulk_Geometry(int *pt, Parameters *params);
  float latice_Potential(const nni ni[7]);

 private:
  int *set_point_type_normals(int *pt, Parameters *params);
};

class Slab_Geometry : public Geometry {
 public:
  Slab_Geometry(int *pt, Parameters *params);
  float latice_Potential(const nni ni[7]);

 private:
  int *set_point_type_normals(int *pt, Parameters *params);
};

class Sphere_Geometry : public Geometry {
 public:
  Sphere_Geometry(int *pt, Parameters *params);
  float latice_Potential(const nni ni[7]);

 private:
  int *set_point_type_normals(int *pt, Parameters *params);
};

class Custom_Geometry : public Geometry {
 public:
  Custom_Geometry(int *pt, Parameters *params);
  float latice_Potential(const nni ni[7]);

 private:
  //~     float *ns;
  int *set_point_type_normals(int *pt, Parameters *params);
};

#endif
