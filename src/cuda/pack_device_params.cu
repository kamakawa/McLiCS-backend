#ifdef USE_CUDA

#include <cuda_runtime.h>
#include <cstring>
#include <string>
#include <vector>
#include <stdexcept>
#include <strings.h> // strcasecmp

#include "../../include/parameters.h"
#include "../../include/cuda/params_device.cuh"
#include "../../include/cuda/anchoring_device.cuh"
#include "../../include/cuda/pack_device_params.cuh"

// =====================
// CUDA error checking
// =====================
static inline void cudaCheck(cudaError_t err, const char* what) {
  if (err != cudaSuccess) {
    std::string msg = std::string("CUDA error on ") + what + ": " + cudaGetErrorString(err);
    throw std::runtime_error(msg);
  }
}

// --------- util: parse boundary type string -> enum ----------
static BoundType parseBound(const char* s) {
  if (!s) return BoundType::Free;
  if (strcasecmp(s, "periodic") == 0) return BoundType::Periodic;
  return BoundType::Free;
}

// --------- util: parse potential string -> enum ----------
static BulkPotentialType parseBulkPotential(const char* s) {
  if (!s) return BulkPotentialType::LL;
  if (strcasecmp(s, "ll") == 0 || strcasecmp(s, "lebwohl-lahser") == 0 || strcasecmp(s, "lebwohl-lasher") == 0)
    return BulkPotentialType::LL;
  if (strcasecmp(s, "ghrl") == 0 || strcasecmp(s, "grun-hess") == 0)
    return BulkPotentialType::GHRL;
  if (strcasecmp(s, "pear") == 0)
    return BulkPotentialType::PEAR;
  return BulkPotentialType::LL;
}

// --------- util: parse anchoring type string -> enum ----------
static AnchoringType parseAnchoring(const std::string& s) {
  if (s == "strong") return AnchoringType::Strong;
  if (s == "strong_ghrl") return AnchoringType::Strong_GHRL;
  if (s == "rp") return AnchoringType::RP;
  if (s == "rp_ghrl") return AnchoringType::RP_GHRL;
  if (s == "fg") return AnchoringType::FG;
  if (s == "fg_ghrl") return AnchoringType::FG_GHRL;
  if (s == "homeotropic") return AnchoringType::Homeotropic;
  if (s == "homeotropic_ghrl") return AnchoringType::Homeotropic_GHRL;
  return AnchoringType::Strong;
}

// --------- pack Parameters -> ParamsDevice ----------
ParamsDevice packParamsDevice(const Parameters& p) {
  ParamsDevice d{};
  d.Nx = p.Nx; d.Ny = p.Ny; d.Nz = p.Nz;

  d.bx = parseBound(p.XBoundtype);
  d.by = parseBound(p.YBoundtype);
  d.bz = parseBound(p.ZBoundtype);

  d.bulkType = parseBulkPotential(p.potential);

  d.A  = p.A;
  d.B1 = p.B1;
  d.B2 = p.B2;
  d.C  = p.C;

  // ============================================================
  // ✅ CRÍTICO: T É A TEMPERATURA REAL (igual CPU).
  // NUNCA colocar invT aqui.
  // ============================================================
  d.T = p.T;

  d.neighbourKind = p.neighbourKind;

  d.ghrl_rho    = p.ghrl_rho;
  d.ghrl_lambda = p.ghrl_lambda;
  d.ghrl_mu     = p.ghrl_mu;
  d.ghrl_nu     = p.ghrl_nu;
  d.ghrl_sigma  = p.ghrl_sigma;

  d.rhoScale    = p.rhoScale;
  d.lambdaScale = p.lambdaScale;
  d.muScale     = p.muScale;
  d.nuScale     = p.nuScale;
  d.sigmaScale  = p.sigmaScale;

  d.elecX = p.elecX; d.elecY = p.elecY; d.elecZ = p.elecZ;
  d.elecA = p.elecA; d.elecE = p.elecE;

  return d;
}

// --------- build anchoring params array from maps ----------
std::vector<AnchoringDeviceParams>
buildAnchoringDeviceParams(const Parameters& p, int numSurfaces) {
  std::vector<AnchoringDeviceParams> out;
  if (numSurfaces <= 0) return out;

  out.resize(numSurfaces);

  for (int id = 0; id < numSurfaces; ++id) {
    AnchoringDeviceParams a{};

    auto itType = p.anchoring_type.find(id);
    a.type = (itType != p.anchoring_type.end()) ? parseAnchoring(itType->second)
                                                : AnchoringType::Strong;

    auto itW = p.W.find(id);
    a.W = (itW != p.W.end()) ? itW->second : 0.0f;

    auto itTh = p.theta_s.find(id);
    a.theta_s = (itTh != p.theta_s.end()) ? itTh->second : 0.0f;

    auto itPh = p.phi_s.find(id);
    a.phi_s = (itPh != p.phi_s.end()) ? itPh->second : 0.0f;

    // (seu rescale por neighbourKind pode ficar aqui quando você for usar surface no CUDA)
    out[id] = a;
  }

  return out;
}

// --------- upload helpers ----------
void uploadParamsToDevice(const ParamsDevice& h, ParamsDevice** d_out) {
  cudaCheck(cudaMalloc((void**)d_out, sizeof(ParamsDevice)), "cudaMalloc(ParamsDevice)");
  cudaCheck(cudaMemcpy(*d_out, &h, sizeof(ParamsDevice), cudaMemcpyHostToDevice),
            "cudaMemcpyHtoD(ParamsDevice)");
}

void uploadAnchoringToDevice(const std::vector<AnchoringDeviceParams>& h,
                             AnchoringDeviceParams** d_out) {
  if (h.empty()) { *d_out = nullptr; return; }
  cudaCheck(cudaMalloc((void**)d_out, sizeof(AnchoringDeviceParams) * h.size()),
            "cudaMalloc(AnchoringDeviceParams)");
  cudaCheck(cudaMemcpy(*d_out, h.data(), sizeof(AnchoringDeviceParams) * h.size(),
                       cudaMemcpyHostToDevice),
            "cudaMemcpyHtoD(AnchoringDeviceParams)");
}

#endif // USE_CUDA