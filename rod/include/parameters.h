#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <gsl/gsl_rng.h>
#include <map>
#include <string>

enum parameter_status { set, unset };

struct Parameters {
  
  struct LatticeParams {
    std::string geometry = "bulk";
    int Nx = 16;
    int Ny = 16;
    int Nz = 16;
    std::string XBoundtype = "free";
    std::string YBoundtype = "free";
    std::string ZBoundtype = "free";
    int (*XBound)(int&, int) = nullptr;
    int (*YBound)(int&, int) = nullptr;
    int (*ZBound)(int&, int) = nullptr;
  } lattice; // Acesso via params.lattice.Nx

  struct SurfaceParams {
    std::string bound_file_name = "";
    std::map<int, float> theta_s;
    std::map<int, float> phi_s;
    std::map<int, float> W;
    std::map<int, std::string> anchoring_type;
  } surface; // Acesso via params.surface.W

  struct PotentialParams {
    std::string potential = "ll";
    float A = 1.0;
    float B1 = 0.04;
    float B2 = 0.4;
    float C = 0.3;
    float Ti = 1.2;
    float Tf = 0.0;
    float dT = -0.02;
    float T = 1.0; // T está aqui.
    float ghrl_rho = 0;
    float ghrl_lambda = 0;
    float ghrl_mu = 0;
    float ghrl_nu = 1;
    float ghrl_sigma = 0;
  } potential; // Acesso via params.potential.T

  struct MCParams {
    std::string evol = "thermal";
    int MCS = 100;
    int MCT = 100;
    int fn = 1;
    int first_file = 0;
  } mc; // Acesso via params.mc.MCS

  struct ICParams {
    std::string ic = "random";
    std::string ic_file = "";
    float phi_0 = 0;
    float theta_0 = 0;
    float p0_i = 0;
    float k11 = 1; // k11, k22, k33 estão aqui.
    float k22 = 1;
    float k33 = 1;
    float p0 = 0;
  } ic; // Acesso via params.ic.k11

  struct NeighbourhoodParams {
    int neighbourKind = 1;
    float rhoScale = 1.0;
  	 float lambdaScale = 1.0;
  	 float muScale = 1.0;
  	 float nuScale = 1.0;
  	 float sigmaScale = 1.0;
  	 float neighbourScale = 1.0;
  } neighbourhood; // Acesso via params.neighbourhood.rhoScale

  struct ElectricFieldParams {
  	 float elecX = 0;
  	 float elecY = 0;
  	 float elecZ = 0;
  	 float elecA = 0;
  	 float elecEi = 0;
  	 float elecEf = 0;
  	 float elecdE = 0;
  	 float elecE = 1;
  } electric; // Acesso via params.electric.elecX
};

#endif