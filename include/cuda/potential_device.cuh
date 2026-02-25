#pragma once
#include "params_device.cuh"

#ifdef __CUDACC__
  #define DINL __device__ __forceinline__
#else
  #define DINL inline
#endif

DINL float bulk_energy_LL(const float ni[3], const float nj[3], const ParamsDevice* p) {
  const float nij = ni[0]*nj[0] + ni[1]*nj[1] + ni[2]*nj[2];
  return -p->A * 1.5f * nij * nij;
}

DINL float bulk_energy_PEAR(const float ni[3], const float nj[3], const ParamsDevice* p, const float rij[3]) {
  const float nij = ni[0]*nj[0] + ni[1]*nj[1] + ni[2]*nj[2];
  const float polar_splay =
      (1.0f + 2.0f * nij + nij * nij) *
      ((nj[0] - ni[0]) * rij[0] + (nj[1] - ni[1]) * rij[1] + (nj[2] - ni[2]) * rij[2]);
  return (-p->A * nij * nij - p->B1 * nij - (p->C * 0.25f) * polar_splay);
}

DINL float bulk_energy_GHRL(const float ni[3], const float nj[3], const ParamsDevice* p, const float rij[3], int nk=1) {
  const float el = (nk == 2 ? p->lambdaScale : 1.0f) * p->ghrl_lambda;
  const float em = (nk == 2 ? p->muScale     : 1.0f) * p->ghrl_mu;
  const float en = (nk == 2 ? p->nuScale     : 1.0f) * p->ghrl_nu;
  const float er = (nk == 2 ? p->rhoScale    : 1.0f) * p->ghrl_rho;
  const float es = (nk == 2 ? p->sigmaScale  : 1.0f) * p->ghrl_sigma;

  const float ai  = ni[0]*rij[0] + ni[1]*rij[1] + ni[2]*rij[2];
  const float aj  = nj[0]*rij[0] + nj[1]*rij[1] + nj[2]*rij[2];
  const float nij = ni[0]*nj[0] + ni[1]*nj[1] + ni[2]*nj[2];

  const float pij = 1.5f * nij * nij - 0.5f;

  const float cross =
      (ni[2]*nj[1] - ni[1]*nj[2]) * rij[0] +
      (ni[0]*nj[2] - ni[2]*nj[0]) * rij[1] +
      (ni[1]*nj[0] - ni[0]*nj[1]) * rij[2];

  const float E1 = (1.5f*ai*ai) + (1.5f*aj*aj) - 1.0f;
  const float sgn = (nij > 0.0f ? 1.0f : -1.0f);

  return ((E1 * (er * pij + el) + em * (ai * aj * nij - (1.0f/9.0f))) + en * pij + es * sgn * cross);
}

DINL float electric_energy(const float ni[3], const ParamsDevice* p) {
  if (p->elecA == 0.0f) return 0.0f;
  const float nDotE = p->elecE * (ni[0]*p->elecX + ni[1]*p->elecY + ni[2]*p->elecZ);
  return -p->elecA * (nDotE * nDotE);
}