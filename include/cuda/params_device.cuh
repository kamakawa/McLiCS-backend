#pragma once

#ifdef __CUDACC__
  #define HD __host__ __device__
#else
  #define HD
#endif

#include <cstdint>

// Boundary types para GPU (em vez de function pointer)
enum class BoundType : int32_t { Free = 0, Periodic = 1 };

// Bulk potential type para GPU (em vez de ponteiro de função em runtime)
enum class BulkPotentialType : int32_t { LL = 0, GHRL = 1, PEAR = 2 };

// Evolve type (caso você queira roteamento no futuro)
enum class EvolveType : int32_t { Thermal = 0, Step = 1, Quench = 2, Electric = 3 };

struct ParamsDevice {
  // lattice
  int32_t Nx, Ny, Nz;
  BoundType bx, by, bz;

  // potential selector
  BulkPotentialType bulkType;

  // potential params
  float A, B1, B2, C;

  // temperature / MC
  float T;
  int32_t neighbourKind;

  // GHRL params + scales
  float ghrl_rho, ghrl_lambda, ghrl_mu, ghrl_nu, ghrl_sigma;
  float rhoScale, lambdaScale, muScale, nuScale, sigmaScale;

  // Electric field
  float elecX, elecY, elecZ;
  float elecA, elecE;
};

// -------- device boundary helpers --------
HD inline int boundary_apply(BoundType bt, int &ii, int NN) {
  if (bt == BoundType::Periodic) {
    ii = (ii + NN) % NN;
    return 1;
  }
  // Free boundary:
  if (ii == -1 || ii == NN) return 0;
  return 1;
}