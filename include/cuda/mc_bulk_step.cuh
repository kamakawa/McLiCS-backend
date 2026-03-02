#pragma once

#ifdef USE_CUDA

#include <cstdint>
#include "../cuda/params_device.cuh"

extern "C" {

// Inicializa (uma vez) buffers no device e faz:
//  - upload do pt (uma vez)
//  - upload do ni inicial (uma vez)
//  - upload do ParamsDevice (uma vez)
void mc_bulk_step_init(
    const float* h_ni,
    const int*   h_pt,
    const ParamsDevice& hp
);

// Atualiza SOMENTE ParamsDevice no device (quando T/campo/etc mudam)
void mc_bulk_step_update_params(const ParamsDevice& hp);

// Executa UM passo Monte Carlo (8 ticks) usando o ni no device,
// e faz apenas o download do ni atualizado para o host (porque seu CPU usa logo após).
void mc_bulk_step_run_download(
    float* h_ni,                 // recebe o ni atualizado (D2H)
    float ang_var,
    std::uint64_t seed,
    int& acceptance_out
);

// (Opcional) Se em algum momento o host modificar ni e você quiser “forçar” o device a refletir:
void mc_bulk_step_upload_ni(const float* h_ni);

// Libera buffers persistentes (opcional; no fim do programa)
void mc_bulk_step_release();

} // extern "C"

#endif // USE_CUDA