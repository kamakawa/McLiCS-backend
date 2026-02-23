#pragma once

#ifdef __CUDACC__
  #define HD __host__ __device__
#else
  #define HD
#endif

#include <cstdint>

enum class AnchoringType : int32_t {
  Strong = 0,
  Strong_GHRL = 1,
  RP = 2,
  RP_GHRL = 3,
  FG = 4,
  FG_GHRL = 5,
  Homeotropic = 6,
  Homeotropic_GHRL = 7
};

struct AnchoringDeviceParams {
  AnchoringType type;
  float W;
  float theta_s; // graus
  float phi_s;   // graus
  // (se no futuro precisar de mais coisas, adiciona aqui)
};