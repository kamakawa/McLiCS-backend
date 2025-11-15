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
#include "../include/anchoring.h"

int Free_Boundary(int &ii, int NN);
int Periodic_Boundary(int &ii, int NN);

class Geometry {
 public:
  Geometry(Parameters *params);
  virtual ~Geometry() = default;

  // Interface pública
  float newman_neighbours(const nni fullni[]);
  float second_neighbours(const nni fullni[]);
  float third_neighbours(const nni fullni[]);
  void Boundary_Init(Parameters *params);
  
  // Getters
  int getNx() const { return Nx; }
  int getNy() const { return Ny; }
  int getNz() const { return Nz; }
  int getNSurfaces() const { return nSurfaces; }
  Parameters* getParams() const { return params; }
  const std::vector<std::unique_ptr<Anchoring>>& getSurfaces() const { return surfaces; }
  
  // Interface para bulk potential
  float computeBulkPotential(float ni[3], float nj[3], Parameters& params, float rij[3], int nk);
  void setBulkPotential(float (*potential_func)(float[3], float[3], Parameters&, float[3], int));

  // Interface abstrata
  virtual int* set_point_type_normals(int* pt, Parameters* params) = 0;
  virtual float lattice_Potential(const nni ni[7]) = 0;
  virtual float Electric_Potential(float ni[3], Parameters& params) = 0;

 protected:
  // Métodos auxiliares protegidos para subclasses
  void addSurface(std::unique_ptr<Anchoring> surface);
  void initializeNS();

 private:
  int Nx, Ny, Nz;
  std::vector<std::unique_ptr<Anchoring>> surfaces;
  int nSurfaces;
  float (*bulk_potential)(float ni[3], float nj[3], Parameters& params, float rij[3], int nk);
  std::unique_ptr<float[]> ns;
  Parameters* params;
};

class Bulk_Geometry : public Geometry {
 public:
  Bulk_Geometry(int* pt, Parameters* params);
  
 private:
  int* set_point_type_normals(int* pt, Parameters* params) override;
  float lattice_Potential(const nni ni[7]) override;
  float Electric_Potential(float ni[3], Parameters& params) override;
};

class Slab_Geometry : public Geometry {
 public:
  Slab_Geometry(int* pt, Parameters* params);
  
 private:
  int* set_point_type_normals(int* pt, Parameters* params) override;
  float lattice_Potential(const nni ni[7]) override;
  float Electric_Potential(float ni[3], Parameters& params) override;
};

class Sphere_Geometry : public Geometry {
 public:
  Sphere_Geometry(int* pt, Parameters* params);
  
 private:
  int* set_point_type_normals(int* pt, Parameters* params) override;
  float lattice_Potential(const nni ni[7]) override;
  float Electric_Potential(float ni[3], Parameters& params) override;
};

class Custom_Geometry : public Geometry {
 public:
  Custom_Geometry(int* pt, Parameters* params);
  
 private:
  int* set_point_type_normals(int* pt, Parameters* params) override;
  float lattice_Potential(const nni ni[7]) override;
  float Electric_Potential(float ni[3], Parameters& params) override;
};

#endif