#include "../include/potential.h"

#include <math.h>
#include <memory>

#include <iostream>

#include "../include/define.h"
#include "../include/parameters.h"

namespace Potential {

float Bulk_Energy_Selinger_Pear(float ni[3], float nj[3], Parameters& params, float rij[3], int nk) {
  float nij, polar_splay;
  nij = ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2];
  polar_splay = +(1 + 2 * nij + nij * nij) * ((nj[0] - ni[0]) * rij[0] + (nj[1] - ni[1]) * rij[1] + (nj[2] - ni[2]) * rij[2]);
  return (-params.A * nij * nij - params.B1 * nij - (params.C * 0.25) * polar_splay);
}

float Bulk_Energy_Lebwohl_Lasher(float ni[3], float nj[3], Parameters& params, float rij[3], int nk) {
  float nij = ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2];
  return -params.A * 1.5 * nij * nij;
}

float Bulk_Energy_GHRL(float ni[3], float nj[3], Parameters& params, float rij[3], int nk) {
  const float el = (nk == 2 ? params.lambdaScale : 1) * params.ghrl_lambda;
  const float em = (nk == 2 ? params.muScale : 1) * params.ghrl_mu;
  const float en = (nk == 2 ? params.nuScale : 1) * params.ghrl_nu;
  const float er = (nk == 2 ? params.rhoScale : 1) * params.ghrl_rho;
  const float es = (nk == 2 ? params.sigmaScale : 1) * params.ghrl_sigma;
  
  float ai = ni[0] * rij[0] + ni[1] * rij[1] + ni[2] * rij[2];
  float aj = nj[0] * rij[0] + nj[1] * rij[1] + nj[2] * rij[2];
  float nij = ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2];
  float pij = 1.5f * nij * nij - 0.5f;
  float cross = (ni[2] * nj[1] - ni[1] * nj[2]) * rij[0] 
              + (ni[0] * nj[2] - ni[2] * nj[0]) * rij[1] 
              + (ni[1] * nj[0] - ni[0] * nj[1]) * rij[2];

  float E1 = (1.5f * ai * ai) + (1.5f * aj * aj) - 1;
  return ((E1 * (er * pij + el) + em * (ai * aj * nij - 1.0f / 9)) + en * pij + es * (nij > 0 ? 1 : -1) * cross);
}

float Electric_Potential(float ni[3], Parameters& params) {
  float nDotE = params.elecE * (ni[0] * params.elecX + ni[1] * params.elecY + ni[2] * params.elecZ);
  return -params.elecA * (nDotE * nDotE);
}

float Bulk_Energy_Selinger_BC(float nix, float niy, float niz, float bix, float biy, float biz, std::unique_ptr<float[]>& ni, std::unique_ptr<float[]>& bi, int i, int j, int k, Parameters& params, std::unique_ptr<int[]>& pt, int[3], int nk) {
  return 0.0f;
}

float nPotential(float nix, float niy, float niz, std::unique_ptr<float[]>& ni, int i, int j, int k, Parameters& params, std::unique_ptr<int[]>& pt, float bulk_pot(float nix, float niy, float niz, std::unique_ptr<float[]>& ni, int i, int j, int k, Parameters& params, std::unique_ptr<int[]>& pt, float rij[3])) {
  return 0.0f;
}

float nbcPotential(float nix, float niy, float niz, float bix, float biy, float biz, std::unique_ptr<float[]>& ni, std::unique_ptr<float[]>& bi, int i, int j, int k, Parameters& params, std::unique_ptr<int[]>& pt,
                   float bulk_pot(float, float, float, float, float, float, std::unique_ptr<float[]>&, std::unique_ptr<float[]>&, int, int, int, Parameters&, std::unique_ptr<int[]>&, int[3])) {
  return 0.0f;
}

} // namespace Potential

int Periodic_Boundary(int &ii, int NN) {
  ii = ((ii + NN) % NN);
  return 1;
}

int Free_Boundary(int &ii, int NN) {
  if (ii == -1 || ii == NN)
    return 0;
  else
    return 1;
}