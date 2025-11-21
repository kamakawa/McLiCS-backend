#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <gsl/gsl_rng.h>
#include <iostream>
#include <vector>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

// Classe base abstrata para definir a geometria do sistema
class Geometry {
 public:
  // Construtor e Destrutor Virtual
  Geometry(Parameters *params);
  virtual ~Geometry() = default;

  // Variaveis de Estado e Configuracao
  int Nx, Ny, Nz;
  int nSurfaces;
  std::vector<class Anchoring *> surfaces;
  float *ns;
  Parameters *params;

  float (*bulk_potential)(float ni[3], float nj[3], Parameters *params, float rij[3], int nk);

  float newman_neighbours(const nni fullni[]);
  float second_nerghbours(const nni fullni[]);
  float third_nerghbours(const nni fullni[]);
  
  void Boundary_Init(Parameters *params);

  virtual int *set_point_type_normals(int *pt, Parameters *params) = 0;
  virtual float latice_Potential(const nni ni[7]) = 0;
};

// Geometria Bulk (Volume livre)
class Bulk_Geometry : public Geometry {
 public:
  Bulk_Geometry(int *pt, Parameters *params);
  
  float latice_Potential(const nni ni[7]) override;

 private:
  int *set_point_type_normals(int *pt, Parameters *params) override;
};

// Geometria do tipo Slab (Placa plana)
class Slab_Geometry : public Geometry {
 public:
  Slab_Geometry(int *pt, Parameters *params);
  
  float latice_Potential(const nni ni[7]) override;

 private:
  int *set_point_type_normals(int *pt, Parameters *params) override;
};

// Geometria Esferica (Gotas ou particulas)
class Sphere_Geometry : public Geometry {
 public:
  Sphere_Geometry(int *pt, Parameters *params);
  
  float latice_Potential(const nni ni[7]) override;

 private:
  int *set_point_type_normals(int *pt, Parameters *params) override;
};

// Geometria Customizada (Definida pelo usuario/arquivo)
class Custom_Geometry : public Geometry {
 public:
  Custom_Geometry(int *pt, Parameters *params);
  
  float latice_Potential(const nni ni[7]) override;

 private:
  int *set_point_type_normals(int *pt, Parameters *params) override;
};

#endif