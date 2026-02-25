#pragma once
#include <cstdint>
#include "params_device.cuh"

#ifdef __cplusplus
extern "C" {
#endif

// Step Monte Carlo para geometria Bulk usando a mesma "arquitetura" do CPU:
// - ParamsDevice (packed do Parameters)
// - bulkType seleciona LL/GHRL/PEAR
// - neighbourKind 1/2/3 (1o/2o/3o vizinhos)
// - inclui electric_energy se elecA!=0
void mc_bulk_step_inplace(
    float* h_ni,            // host: 3*N floats
    int*   h_pt,            // host: pelo menos N ints (pt[idx])
    const ParamsDevice& hp, // host params packed
    float ang_var,
    std::uint64_t seed,
    int& acceptance_out
);

#ifdef __cplusplus
}
#endif