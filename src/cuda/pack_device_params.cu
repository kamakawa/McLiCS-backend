#include <cuda_runtime.h>
#include <cstring>
#include <string>
#include <vector>
#include <stdexcept>

#include "../include/parameters.h"
#include "../include/cuda/params_device.cuh"
#include "../include/cuda/anchoring_device.cuh"

// --------- util: parse boundary type string -> enum ----------
static BoundType parseBound(const char* s) {
  if (!s) return BoundType::Free;
  if (strcasecmp(s, "periodic") == 0) return BoundType::Periodic;
  return BoundType::Free; // default compatível com seu Parameters
}

// --------- util: parse potential string -> enum ----------
static BulkPotentialType parseBulkPotential(const char* s) {
  if (!s) return BulkPotentialType::LL;
  if (strcasecmp(s, "ll") == 0 || strcasecmp(s, "lebwohl-lahser") == 0) return BulkPotentialType::LL;
  if (strcasecmp(s, "ghrl") == 0 || strcasecmp(s, "grun-hess") == 0) return BulkPotentialType::GHRL;
  if (strcasecmp(s, "pear") == 0) return BulkPotentialType::PEAR;
  return BulkPotentialType::LL; // fallback conservador
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
  // default seguro:
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

  d.A = p.A; d.B1 = p.B1; d.B2 = p.B2; d.C = p.C;

  d.T = p.T;
  d.neighbourKind = p.neighbourKind;

  d.ghrl_rho = p.ghrl_rho;
  d.ghrl_lambda = p.ghrl_lambda;
  d.ghrl_mu = p.ghrl_mu;
  d.ghrl_nu = p.ghrl_nu;
  d.ghrl_sigma = p.ghrl_sigma;

  d.rhoScale = p.rhoScale;
  d.lambdaScale = p.lambdaScale;
  d.muScale = p.muScale;
  d.nuScale = p.nuScale;
  d.sigmaScale = p.sigmaScale;

  d.elecX = p.elecX; d.elecY = p.elecY; d.elecZ = p.elecZ;
  d.elecA = p.elecA; d.elecE = p.elecE;

  return d;
}

// --------- build anchoring params array from maps ----------
std::vector<AnchoringDeviceParams> buildAnchoringDeviceParams(const Parameters& p, int numSurfaces) {
  std::vector<AnchoringDeviceParams> out;
  out.resize(numSurfaces);

  for (int id = 0; id < numSurfaces; ++id) {
    AnchoringDeviceParams a{};

    // type
    auto itType = p.anchoring_type.find(id);
    if (itType != p.anchoring_type.end()) a.type = parseAnchoring(itType->second);
    else a.type = AnchoringType::Strong;

    // W
    auto itW = p.W.find(id);
    a.W = (itW != p.W.end()) ? itW->second : 0.0f;

    // theta/phi (em graus, como seu código usa)
    auto itTh = p.theta_s.find(id);
    a.theta_s = (itTh != p.theta_s.end()) ? itTh->second : 0.0f;

    auto itPh = p.phi_s.find(id);
    a.phi_s = (itPh != p.phi_s.end()) ? itPh->second : 0.0f;

    out[id] = a;
  }
  return out;
}

// --------- upload helpers ----------
void uploadParamsToDevice(const ParamsDevice& h, ParamsDevice** d_out) {
  cudaMalloc((void**)d_out, sizeof(ParamsDevice));
  cudaMemcpy(*d_out, &h, sizeof(ParamsDevice), cudaMemcpyHostToDevice);
}

void uploadAnchoringToDevice(const std::vector<AnchoringDeviceParams>& h, AnchoringDeviceParams** d_out) {
  if (h.empty()) { *d_out = nullptr; return; }
  cudaMalloc((void**)d_out, sizeof(AnchoringDeviceParams) * h.size());
  cudaMemcpy(*d_out, h.data(), sizeof(AnchoringDeviceParams) * h.size(), cudaMemcpyHostToDevice);
}