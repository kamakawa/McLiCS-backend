#pragma once
#ifdef USE_CUDA

#include <vector>
#include "../parameters.h"
#include "../cuda/params_device.cuh"
#include "../cuda/anchoring_device.cuh"

// Converte Parameters (CPU) -> ParamsDevice (GPU)
ParamsDevice packParamsDevice(const Parameters& p);

// Constrói array de parâmetros de ancoramento (se/quando usar)
std::vector<AnchoringDeviceParams>
buildAnchoringDeviceParams(const Parameters& p, int numSurfaces);

// Upload helpers
void uploadParamsToDevice(const ParamsDevice& h, ParamsDevice** d_out);
void uploadAnchoringToDevice(const std::vector<AnchoringDeviceParams>& h,
                             AnchoringDeviceParams** d_out);

#endif // USE_CUDA