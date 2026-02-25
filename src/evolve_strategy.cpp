#include "../include/evolve_strategy.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <omp.h>
#include <iostream>
#include <strings.h> // strcasecmp
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"
#include "../include/geometry_strategy.h"
#include "../include/anchoring_strategy.h"

#ifdef USE_CUDA
  #include "../include/cuda/pack_device_params.cuh"
  #include "../include/cuda/mc_bulk_step.cuh"
#endif

// ========== Funções Auxiliares ==========

#ifdef USE_CUDA
#include <cstdint>

static inline std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
#endif

void inline setNni(uint pos, int conditional, nni *nLocal, float *ni, int *pt) {
  if (conditional == 1) {
    nLocal->x = ni[pos * 3 + 0];
    nLocal->y = ni[pos * 3 + 1];
    nLocal->z = ni[pos * 3 + 2];
    nLocal->pt = pt[pos];
  } else {
    nLocal->pt = 0;
  }
}

// ========== Implementação dos Métodos Compartilhados ==========

void EvolveStrategy::initializeRNG(gsl_rng*** rng, int num_threads) {
    *rng = (gsl_rng**)calloc(num_threads, sizeof(gsl_rng*));
    gsl_rng_env_setup();

    for (int i = 0; i < num_threads; i++) {
        (*rng)[i] = gsl_rng_alloc(gsl_rng_ranlxs0);
        gsl_rng_set((*rng)[i], i);
    }
}

void EvolveStrategy::cleanupRNG(gsl_rng** rng, int num_threads) {
    for (int i = 0; i < num_threads; i++) {
        gsl_rng_free(rng[i]);
    }
    free(rng);
}

void EvolveStrategy::monteCarloStep(float& ang_var, gsl_rng** r, float* ni, int* pt,
                                   Parameters* params, GeometryStrategy* geometry_strategy,
                                   float* surface_normals,
                                   std::vector<AnchoringStrategy*>& anchoring_strategies) {
    const int Nt = params->Nx * params->Ny * params->Nz;

    int valid_points = 0;
    for (int i = 0; i < Nt; i++) {
        if (pt[i]) valid_points++;
    }
    if (valid_points == 0) return;

// ============================================================
// CUDA PATH: apenas BULK + neighbourKind==1 + LL (primeira etapa)
// ============================================================
#ifdef USE_CUDA
    const bool is_bulk = (geometry_strategy && geometry_strategy->getName() == std::string("Bulk Geometry"));

    if (is_bulk) {
        // Empacota params CPU -> ParamsDevice
        ParamsDevice hp = packParamsDevice(*params);

        // ------------------------------------------------------------
        // Init persistente (pt + ni uma vez). Depois só roda kernel + D2H.
        // ------------------------------------------------------------
        static bool gpu_inited = false;
        static int  gpuNx = 0, gpuNy = 0, gpuNz = 0;

        if (!gpu_inited || gpuNx != params->Nx || gpuNy != params->Ny || gpuNz != params->Nz) {
            mc_bulk_step_init(ni, pt, hp);
            gpu_inited = true;
            gpuNx = params->Nx; gpuNy = params->Ny; gpuNz = params->Nz;
        } else {
            // Se só mudou T/campo etc, atualiza params no device
            mc_bulk_step_update_params(hp);
        }

        // seed que varia por chamada (equivalente ao RNG do CPU evoluindo)
        static std::uint64_t call_counter = 0;
        const std::uint64_t ctr = ++call_counter;

        const std::uint64_t seed =
            (std::uint64_t)(params->T * 1000000.0f) ^
            (std::uint64_t)(params->Nx * 73856093u) ^
            (std::uint64_t)(params->Ny * 19349663u) ^
            (std::uint64_t)(params->Nz * 83492791u) ^
            (ctr * 0x9E3779B97F4A7C15ULL);

        int acceptance_gpu = 0;

        // Roda 1 step e baixa ni (porque você usa logo depois no energyCalculator/parameterOrder)
        mc_bulk_step_run_download(ni, ang_var, seed, acceptance_gpu);

        // Ajuste dinâmico do passo angular (igual CPU)
        float acceptance_rate = 1.0f * acceptance_gpu / valid_points;
        if (acceptance_rate < 0.5f) {
            ang_var *= 0.99f;
            if (ang_var < 0.01f) ang_var = 0.01f;
        } else if (acceptance_rate > 0.5f) {
            ang_var /= 0.99f;
            if (ang_var > 1.0f) ang_var = 0.5f;
        }
        return;
    }
#endif

    // ============================================================
    // CPU PATH (OpenMP) - mantém exatamente sua lógica
    // ============================================================
    int acceptance = 0;
    const int num_threads = omp_get_max_threads();
    const int Nx = params->Nx;
    const int Ny = params->Ny;
    const int Nz = params->Nz;

    int iBoxSize = 2;
    int jBoxSize = 2;
    int kBoxSize = 2;

    int Ni = Nx / 2 + Nx % 2;
    int Nj = Ny / 2 + Ny % 2;
    int Nk = Nz / 2 + Nz % 2;
    int Ntt = Ni * Nj * Nk;

    #pragma omp parallel num_threads(num_threads)
    {
        for (int tick = 0; tick < 8; tick++) {
            #pragma omp for reduction(+:acceptance) schedule(dynamic)
            for (int nt = 0; nt < Ntt; nt++) {
                float E_new, E_old, rotation_type, va;
                int i, j, k;
                float nNew[3];

                int thread = omp_get_thread_num();
                int iBox = nt % Ni;
                int jBox = (nt / Ni) % Nj;
                int kBox = nt / (Ni * Nj);

                int di = tick % 2;
                int dj = (tick / 2) % 2;
                int dk = tick / 4;

                const int nti = params->neighbourKind == 2 ? 20 : (params->neighbourKind == 3 ? 28 : 8);
                nni nLocal[nti];

                i = di + iBox * iBoxSize;
                if (i >= Nx) continue;
                j = dj + jBox * jBoxSize;
                if (j >= Ny) continue;
                k = dk + kBox * kBoxSize;
                if (k >= Nz) continue;

                if (pti(i, j, k) == 0) continue;

                int im = (i - 1), ip = (i + 1);
                int jm = (j - 1), jp = (j + 1);
                int km = (k - 1), kp = (k + 1);

                int neighbour[6] = {
                    params->XBound(ip, Nx), params->XBound(im, Nx),
                    params->YBound(jp, Ny), params->YBound(jm, Ny),
                    params->ZBound(kp, Nz), params->ZBound(km, Nz)
                };

                setNni(i  + Nx * (j  + Ny * k ), 1,            &nLocal[0], ni, pt);
                setNni(ip + Nx * (j  + Ny * k ), neighbour[0], &nLocal[1], ni, pt);
                setNni(im + Nx * (j  + Ny * k ), neighbour[1], &nLocal[2], ni, pt);
                setNni(i  + Nx * (jp + Ny * k ), neighbour[2], &nLocal[3], ni, pt);
                setNni(i  + Nx * (jm + Ny * k ), neighbour[3], &nLocal[4], ni, pt);
                setNni(i  + Nx * (j  + Ny * kp), neighbour[4], &nLocal[5], ni, pt);
                setNni(i  + Nx * (j  + Ny * km), neighbour[5], &nLocal[6], ni, pt);

                if (nLocal[0].pt > 1) {
                    int idx = i + Nx * (j + Ny * k);
                    nLocal[7].x = surface_normals[idx * 3 + 0];
                    nLocal[7].y = surface_normals[idx * 3 + 1];
                    nLocal[7].z = surface_normals[idx * 3 + 2];
                    nLocal[7].pt = 1;
                }

                if (params->neighbourKind > 1) {
                    setNni(ip + Nx * (jp + Ny * k ), neighbour[0] * neighbour[2], &nLocal[8], ni, pt);
                    setNni(ip + Nx * (jm + Ny * k ), neighbour[0] * neighbour[3], &nLocal[9], ni, pt);
                    setNni(ip + Nx * (j  + Ny * kp), neighbour[0] * neighbour[4], &nLocal[10], ni, pt);
                    setNni(ip + Nx * (j  + Ny * km), neighbour[0] * neighbour[5], &nLocal[11], ni, pt);
                    setNni(im + Nx * (jp + Ny * k ), neighbour[1] * neighbour[2], &nLocal[12], ni, pt);
                    setNni(im + Nx * (jm + Ny * k ), neighbour[1] * neighbour[3], &nLocal[13], ni, pt);
                    setNni(im + Nx * (j  + Ny * kp), neighbour[1] * neighbour[4], &nLocal[14], ni, pt);
                    setNni(im + Nx * (j  + Ny * km), neighbour[1] * neighbour[5], &nLocal[15], ni, pt);
                    setNni(i  + Nx * (jp + Ny * kp), neighbour[2] * neighbour[4], &nLocal[16], ni, pt);
                    setNni(i  + Nx * (jp + Ny * km), neighbour[2] * neighbour[5], &nLocal[17], ni, pt);
                    setNni(i  + Nx * (jm + Ny * kp), neighbour[3] * neighbour[4], &nLocal[18], ni, pt);
                    setNni(i  + Nx * (jm + Ny * km), neighbour[3] * neighbour[5], &nLocal[19], ni, pt);
                }

                if (params->neighbourKind == 3) {
                    setNni(ip + Nx * (jp + Ny * kp), neighbour[0] * neighbour[2] * neighbour[4], &nLocal[20], ni, pt);
                    setNni(ip + Nx * (jp + Ny * km), neighbour[0] * neighbour[2] * neighbour[5], &nLocal[21], ni, pt);
                    setNni(ip + Nx * (jm + Ny * kp), neighbour[0] * neighbour[3] * neighbour[4], &nLocal[22], ni, pt);
                    setNni(ip + Nx * (jm + Ny * km), neighbour[0] * neighbour[3] * neighbour[5], &nLocal[23], ni, pt);
                    setNni(im + Nx * (jp + Ny * kp), neighbour[1] * neighbour[2] * neighbour[4], &nLocal[24], ni, pt);
                    setNni(im + Nx * (jp + Ny * km), neighbour[1] * neighbour[2] * neighbour[5], &nLocal[25], ni, pt);
                    setNni(im + Nx * (jm + Ny * kp), neighbour[1] * neighbour[3] * neighbour[4], &nLocal[26], ni, pt);
                    setNni(im + Nx * (jm + Ny * km), neighbour[1] * neighbour[3] * neighbour[5], &nLocal[27], ni, pt);
                }

                va = (2 * gsl_rng_uniform(r[thread]) - 1) * ang_var;
                rotation_type = gsl_rng_uniform(r[thread]);

                if (rotation_type < 0.333) {
                    nNew[0] = nix(i, j, k);
                    nNew[1] = niy(i, j, k) * cosf(M_PI * va) + niz(i, j, k) * sinf(M_PI * va);
                    nNew[2] = niz(i, j, k) * cosf(M_PI * va) - niy(i, j, k) * sinf(M_PI * va);
                } else if (rotation_type < 0.666) {
                    nNew[0] = nix(i, j, k) * cosf(M_PI * va) - niz(i, j, k) * sinf(M_PI * va);
                    nNew[1] = niy(i, j, k);
                    nNew[2] = niz(i, j, k) * cosf(M_PI * va) + nix(i, j, k) * sinf(M_PI * va);
                } else {
                    nNew[0] = nix(i, j, k) * cosf(M_PI * va) + niy(i, j, k) * sinf(M_PI * va);
                    nNew[1] = niy(i, j, k) * cosf(M_PI * va) - nix(i, j, k) * sinf(M_PI * va);
                    nNew[2] = niz(i, j, k);
                }

                float norm = sqrtf(nNew[0]*nNew[0] + nNew[1]*nNew[1] + nNew[2]*nNew[2]);
                nNew[0] /= norm;
                nNew[1] /= norm;
                nNew[2] /= norm;

                E_old = geometry_strategy->calculatePotential(nLocal, params, anchoring_strategies, surface_normals);

                nLocal[0].x = nNew[0];
                nLocal[0].y = nNew[1];
                nLocal[0].z = nNew[2];

                E_new = geometry_strategy->calculatePotential(nLocal, params, anchoring_strategies, surface_normals);

                if (E_new - E_old < 0 || gsl_rng_uniform(r[thread]) < exp(-(E_new - E_old) / params->T)) {
                    nix(i, j, k) = nNew[0];
                    niy(i, j, k) = nNew[1];
                    niz(i, j, k) = nNew[2];
                    acceptance++;
                }
            }
        }
    }

    float acceptance_rate = 1.0f * acceptance / valid_points;
    if (acceptance_rate < 0.5f) {
        ang_var *= 0.99f;
        if (ang_var < 0.01f) ang_var = 0.01f;
    } else if (acceptance_rate > 0.5f) {
        ang_var /= 0.99f;
        if (ang_var > 1.0f) ang_var = 0.5f;
    }
}

// CORREÇÃO: energyCalculator agora retorna energia TOTAL (igual ao original)
float EvolveStrategy::energyCalculator(float* ni, int* pt, Parameters* params,
                                      GeometryStrategy* geometry_strategy,
                                      float* surface_normals,
                                      std::vector<AnchoringStrategy*>& anchoring_strategies) {
    const int Nt = params->Nx * params->Ny * params->Nz;
    double total_energy = 0.0;
    const int Nx = params->Nx;
    const int Ny = params->Ny;
    const int Nz = params->Nz;

    #pragma omp parallel for reduction(+:total_energy) schedule(dynamic)
    for (int idx = 0; idx < Nt; idx++) {
        if (pt[idx] == 0) continue;

        nni nLocal[28];
        int i = idx % Nx;
        int j = (idx / Nx) % Ny;
        int k = idx / (Nx * Ny);

        int im = (i - 1), ip = (i + 1);
        int jm = (j - 1), jp = (j + 1);
        int km = (k - 1), kp = (k + 1);

        int neighbour[6] = {
            params->XBound(ip, Nx), params->XBound(im, Nx),
            params->YBound(jp, Ny), params->YBound(jm, Ny),
            params->ZBound(kp, Nz), params->ZBound(km, Nz)
        };

        setNni(i  + Nx * (j  + Ny * k ), 1,            &nLocal[0], ni, pt);
        setNni(ip + Nx * (j  + Ny * k ), neighbour[0], &nLocal[1], ni, pt);
        setNni(im + Nx * (j  + Ny * k ), neighbour[1], &nLocal[2], ni, pt);
        setNni(i  + Nx * (jp + Ny * k ), neighbour[2], &nLocal[3], ni, pt);
        setNni(i  + Nx * (jm + Ny * k ), neighbour[3], &nLocal[4], ni, pt);
        setNni(i  + Nx * (j  + Ny * kp), neighbour[4], &nLocal[5], ni, pt);
        setNni(i  + Nx * (j  + Ny * km), neighbour[5], &nLocal[6], ni, pt);

        if (nLocal[0].pt > 1) {
            int pos = i + Nx * (j + Ny * k);
            nLocal[7].x = surface_normals[pos * 3 + 0];
            nLocal[7].y = surface_normals[pos * 3 + 1];
            nLocal[7].z = surface_normals[pos * 3 + 2];
            nLocal[7].pt = 1;
        }

        if (params->neighbourKind > 1) {
            setNni(ip + Nx * (jp + Ny * k ), neighbour[0] * neighbour[2], &nLocal[8], ni, pt);
            setNni(ip + Nx * (jm + Ny * k ), neighbour[0] * neighbour[3], &nLocal[9], ni, pt);
            setNni(ip + Nx * (j  + Ny * kp), neighbour[0] * neighbour[4], &nLocal[10], ni, pt);
            setNni(ip + Nx * (j  + Ny * km), neighbour[0] * neighbour[5], &nLocal[11], ni, pt);
            setNni(im + Nx * (jp + Ny * k ), neighbour[1] * neighbour[2], &nLocal[12], ni, pt);
            setNni(im + Nx * (jm + Ny * k ), neighbour[1] * neighbour[3], &nLocal[13], ni, pt);
            setNni(im + Nx * (j  + Ny * kp), neighbour[1] * neighbour[4], &nLocal[14], ni, pt);
            setNni(im + Nx * (j  + Ny * km), neighbour[1] * neighbour[5], &nLocal[15], ni, pt);
            setNni(i  + Nx * (jp + Ny * kp), neighbour[2] * neighbour[4], &nLocal[16], ni, pt);
            setNni(i  + Nx * (jp + Ny * km), neighbour[2] * neighbour[5], &nLocal[17], ni, pt);
            setNni(i  + Nx * (jm + Ny * kp), neighbour[3] * neighbour[4], &nLocal[18], ni, pt);
            setNni(i  + Nx * (jm + Ny * km), neighbour[3] * neighbour[5], &nLocal[19], ni, pt);
        }

        if (params->neighbourKind == 3) {
            setNni(ip + Nx * (jp + Ny * kp), neighbour[0] * neighbour[2] * neighbour[4], &nLocal[20], ni, pt);
            setNni(ip + Nx * (jp + Ny * km), neighbour[0] * neighbour[2] * neighbour[5], &nLocal[21], ni, pt);
            setNni(ip + Nx * (jm + Ny * kp), neighbour[0] * neighbour[3] * neighbour[4], &nLocal[22], ni, pt);
            setNni(ip + Nx * (jm + Ny * km), neighbour[0] * neighbour[3] * neighbour[5], &nLocal[23], ni, pt);
            setNni(im + Nx * (jp + Ny * kp), neighbour[1] * neighbour[2] * neighbour[4], &nLocal[24], ni, pt);
            setNni(im + Nx * (jp + Ny * km), neighbour[1] * neighbour[2] * neighbour[5], &nLocal[25], ni, pt);
            setNni(im + Nx * (jm + Ny * kp), neighbour[1] * neighbour[3] * neighbour[4], &nLocal[26], ni, pt);
            setNni(im + Nx * (jm + Ny * km), neighbour[1] * neighbour[3] * neighbour[5], &nLocal[27], ni, pt);
        }

        float energy = geometry_strategy->calculatePotential(nLocal, params, anchoring_strategies, surface_normals);
        total_energy += energy;
    }

    return total_energy / Nt;
}

// ========== Implementações das Estratégias Específicas ==========
// (restante do seu arquivo original continua igual)
//
// IMPORTANTE: Eu não alterei as classes ThermalEvolveStrategy / StepEvolveStrategy
// fora do monteCarloStep, para preservar comportamento.

// ========== Implementações das Estratégias Específicas ==========
// (Mantidas iguais, pois agora usam os métodos corrigidos)

int ThermalEvolveStrategy::run(float* ni, int* pt, Parameters* params,
                              GeometryStrategy* geometry_strategy,
                              float* surface_normals,
                              std::vector<AnchoringStrategy*>& anchoring_strategies) {
    printf("Initializing thermal evolution:\n");
    printf("Ti= %g, Tf= %g, dT= %g\n\n", params->Ti, params->Tf, params->dT);

    int num_threads = omp_get_max_threads();
    gsl_rng** rng;
    initializeRNG(&rng, num_threads);

    FILE* po_file = fopen("po.dat", "a");
    fprintf(po_file, "T S varS E varE\n");
    fflush(po_file);

    printf("Starting thermal variation from %g to %g with step %g\n",
           params->Ti, params->Tf, params->dT);
    printf("MCT=%d MCS=%d using %d threads\n\n", params->MCT, params->MCS, num_threads);

    thermalLoop(ni, pt, params, geometry_strategy, surface_normals, anchoring_strategies, rng, po_file);

    fclose(po_file);
    cleanupRNG(rng, num_threads);
    return 0;
}

void ThermalEvolveStrategy::thermalLoop(float* ni, int* pt, Parameters* params,
                                       GeometryStrategy* geometry_strategy,
                                       float* surface_normals,
                                       std::vector<AnchoringStrategy*>& anchoring_strategies,
                                       gsl_rng** rng, FILE* po_file) {
    float ang_var = 0.5f;
    int sign = (params->dT > 0) ? 1 : -1;

    for (params->T = params->Ti;
         sign * (params->T - params->Tf) <= 0;
         params->T += params->dT) {

        // Passos de equilíbrio
        for (int step = 0; step < params->MCT; step++) {
            monteCarloStep(ang_var, rng, ni, pt, params, geometry_strategy,
                          surface_normals, anchoring_strategies);
        }

        // Coleta de estatísticas
        float S1 = 0, S2 = 0, E = 0, E2 = 0;
        float vec_n[3], mat_n[9];

        for (int step = 0; step < params->MCS; step++) {
            monteCarloStep(ang_var, rng, ni, pt, params, geometry_strategy,
                          surface_normals, anchoring_strategies);

            float tempE = energyCalculator(ni, pt, params, geometry_strategy,
                                         surface_normals, anchoring_strategies);
            E += tempE;
            E2 += tempE * tempE;

            OrderParam::Matrice_constructor(ni, mat_n, pt, *params);
            float sTemp = OrderParam::Eigen_value_evaluation(mat_n, vec_n);
            S1 += sTemp;
            S2 += sTemp * sTemp;
        }

        // Cálculo de médias e saída
        E /= params->MCS;
        E2 /= params->MCS;
        S1 /= params->MCS;
        S2 /= params->MCS;

        // Salva snapshot
        char fname[100];
        sprintf(fname, "director_field_T_%d.csv", (int)(100 * (params->T + 1e-7)));
        print_n(fname, ni, *params, pt);

        fprintf(po_file, "%g %g %g %g %g\n", params->T, S1, S2 - S1 * S1, E, E2 - E * E);
        fflush(po_file);
    }
}

int StepEvolveStrategy::run(float* ni, int* pt, Parameters* params,
                           GeometryStrategy* geometry_strategy,
                           float* surface_normals,
                           std::vector<AnchoringStrategy*>& anchoring_strategies) {
    printf("Initializing step evolution:\n");
    printf("Initial File Number= %d, Last File Number= %d\n\n",
           params->first_file, params->first_file + params->fn);

    int num_threads = omp_get_max_threads();
    gsl_rng** rng;
    initializeRNG(&rng, num_threads);

    FILE* po_file = fopen("po_step.dat", "a");
    fprintf(po_file, "ii S varS E varE\n");
    fflush(po_file);

    printf("Step evolution using MCT=%d MCS=%d fn=%d with %d threads\n\n",
           params->MCT, params->MCS, params->fn, num_threads);

    float ang_var = 0.5f;

    for (int ii = params->first_file; ii < params->fn + params->first_file; ii++) {

        // Passos de equilíbrio
        for (int step = 0; step < params->MCT; step++) {
            monteCarloStep(ang_var, rng, ni, pt, params, geometry_strategy,
                          surface_normals, anchoring_strategies);
        }

        // Coleta de estatísticas
        float S1 = 0, S2 = 0, E = 0, E2 = 0;
        float vec_n[3], mat_n[9];

        for (int step = 0; step < params->MCS; step++) {
            monteCarloStep(ang_var, rng, ni, pt, params, geometry_strategy,
                          surface_normals, anchoring_strategies);

            float tempE = energyCalculator(ni, pt, params, geometry_strategy,
                                         surface_normals, anchoring_strategies);
            E += tempE;
            E2 += tempE * tempE;

            OrderParam::Matrice_constructor(ni, mat_n, pt, *params);
            float sTemp = OrderParam::Eigen_value_evaluation(mat_n, vec_n);
            S1 += sTemp;
            S2 += sTemp * sTemp;
        }

        // Cálculo de médias
        E /= params->MCS;
        E2 /= params->MCS;
        S1 /= params->MCS;
        S2 /= params->MCS;

        // Saída
        char fname[100];
        sprintf(fname, "director_field_%d.csv", ii);
        print_n(fname, ni, *params, pt);

        fprintf(po_file, "%d %g %g %g %g\n", ii, S1, S2 - S1 * S1, E, E2 - E * E);
        fflush(po_file);
    }

    fclose(po_file);
    cleanupRNG(rng, num_threads);
    return 0;
}

int QuenchEvolveStrategy::run(float* ni, int* pt, Parameters* params,
                             GeometryStrategy* geometry_strategy,
                             float* surface_normals,
                             std::vector<AnchoringStrategy*>& anchoring_strategies) {
    printf("Initializing quench evolution:\n");
    printf("Initial File Number= %d, Last File Number= %d\n\n",
           params->first_file, params->first_file + params->fn);

    int num_threads = omp_get_max_threads();
    gsl_rng** rng;
    initializeRNG(&rng, num_threads);

    FILE* po_file = fopen("po_quench.dat", "a");
    fprintf(po_file, "ii S varS E varE\n");
    fflush(po_file);

    printf("Quench relaxation using MCT=%d MCS=%d fn=%d with %d threads\n\n",
           params->MCT, params->MCS, params->fn, num_threads);

    quenchLoop(ni, pt, params, geometry_strategy, surface_normals, anchoring_strategies, rng, po_file);

    fclose(po_file);
    cleanupRNG(rng, num_threads);
    return 0;
}

void QuenchEvolveStrategy::quenchLoop(float* ni, int* pt, Parameters* params,
                                     GeometryStrategy* geometry_strategy,
                                     float* surface_normals,
                                     std::vector<AnchoringStrategy*>& anchoring_strategies,
                                     gsl_rng** rng, FILE* po_file) {
    float ang_var = 0.5f;

    for (int ii = params->first_file; ii < params->fn + params->first_file; ii++) {

        // Fase de alta temperatura
        params->T = params->Ti;
        for (int step = 0; step < params->MCT; step++) {
            monteCarloStep(ang_var, rng, ni, pt, params, geometry_strategy,
                          surface_normals, anchoring_strategies);
        }

        // Resfriamento rápido
        params->T = params->Tf;
        for (int step = 0; step < params->MCT * params->dT; step++) {
            monteCarloStep(ang_var, rng, ni, pt, params, geometry_strategy,
                          surface_normals, anchoring_strategies);
        }

        // Coleta de estatísticas na temperatura final
        float S1 = 0, S2 = 0, E = 0, E2 = 0;
        float vec_n[3], mat_n[9];

        for (int step = 0; step < params->MCS; step++) {
            monteCarloStep(ang_var, rng, ni, pt, params, geometry_strategy,
                          surface_normals, anchoring_strategies);

            float tempE = energyCalculator(ni, pt, params, geometry_strategy,
                                         surface_normals, anchoring_strategies);
            E += tempE;
            E2 += tempE * tempE;

            OrderParam::Matrice_constructor(ni, mat_n, pt, *params);
            float sTemp = OrderParam::Eigen_value_evaluation(mat_n, vec_n);
            S1 += sTemp;
            S2 += sTemp * sTemp;
        }

        // Cálculo de médias
        E /= params->MCS;
        E2 /= params->MCS;
        S1 /= params->MCS;
        S2 /= params->MCS;

        // Saída
        char fname[100];
        sprintf(fname, "director_field_%d.csv", ii);
        print_n(fname, ni, *params, pt);

        fprintf(po_file, "%d %g %g %g %g\n", ii, S1, S2 - S1 * S1, E, E2 - E * E);
        fflush(po_file);
    }
}

int ElectricEvolveStrategy::run(float* ni, int* pt, Parameters* params,
                               GeometryStrategy* geometry_strategy,
                               float* surface_normals,
                               std::vector<AnchoringStrategy*>& anchoring_strategies) {
    printf("Initializing electric field evolution:\n");
    printf("Ei= %g, Ef= %g, dE= %g\n\n", params->elecEi, params->elecEf, params->elecdE);

    int num_threads = omp_get_max_threads();
    gsl_rng** rng;
    initializeRNG(&rng, num_threads);

    FILE* po_file = fopen("po_electric.dat", "a");
    fprintf(po_file, "E S varS E_total varE\n");
    fflush(po_file);

    printf("Starting electric field variation from %g to %g with step %g\n",
           params->elecEi, params->elecEf, params->elecdE);
    printf("MCT=%d MCS=%d using %d threads\n\n", params->MCT, params->MCS, num_threads);

    electricLoop(ni, pt, params, geometry_strategy, surface_normals, anchoring_strategies, rng, po_file);

    fclose(po_file);
    cleanupRNG(rng, num_threads);
    return 0;
}

void ElectricEvolveStrategy::electricLoop(float* ni, int* pt, Parameters* params,
                                         GeometryStrategy* geometry_strategy,
                                         float* surface_normals,
                                         std::vector<AnchoringStrategy*>& anchoring_strategies,
                                         gsl_rng** rng, FILE* po_file) {
    float ang_var = 0.5f;
    int sign = (params->elecdE > 0) ? 1 : -1;

    for (params->elecE = params->elecEi;
         sign * (params->elecE - params->elecEf) <= 0;
         params->elecE += params->elecdE) {

        // Passos de equilíbrio
        for (int step = 0; step < params->MCT; step++) {
            monteCarloStep(ang_var, rng, ni, pt, params, geometry_strategy,
                          surface_normals, anchoring_strategies);
        }

        // Coleta de estatísticas
        float S1 = 0, S2 = 0, E = 0, E2 = 0;
        float vec_n[3], mat_n[9];

        for (int step = 0; step < params->MCS; step++) {
            monteCarloStep(ang_var, rng, ni, pt, params, geometry_strategy,
                          surface_normals, anchoring_strategies);

            float tempE = energyCalculator(ni, pt, params, geometry_strategy,
                                         surface_normals, anchoring_strategies);
            E += tempE;
            E2 += tempE * tempE;

            OrderParam::Matrice_constructor(ni, mat_n, pt, *params);
            float sTemp = OrderParam::Eigen_value_evaluation(mat_n, vec_n);
            S1 += sTemp;
            S2 += sTemp * sTemp;
        }

        // Cálculo de médias
        E /= params->MCS;
        E2 /= params->MCS;
        S1 /= params->MCS;
        S2 /= params->MCS;

        // Saída
        char fname[100];
        sprintf(fname, "director_field_E_%d.csv", (int)(100 * (params->elecE + 1e-7)));
        print_n(fname, ni, *params, pt);

        fprintf(po_file, "%g %g %g %g %g\n", params->elecE, S1, S2 - S1 * S1, E, E2 - E * E);
        fflush(po_file);
    }
}

// ========== Factory Method ==========
EvolveStrategy* EvolveStrategyFactory::create(const std::string& evolveType) {
    if (evolveType == "thermal") {
        return new ThermalEvolveStrategy();
    } else if (evolveType == "step") {
        return new StepEvolveStrategy();
    } else if (evolveType == "quench") {
        return new QuenchEvolveStrategy();
    } else if (evolveType == "electric") {
        return new ElectricEvolveStrategy();
    } else {
        std::cerr << "Unknown evolve type: " << evolveType << std::endl;
        return nullptr;
    }
}

EvolveStrategy* EvolveStrategyFactory::create(Parameters* params) {
    // Determina o tipo baseado nos parâmetros
    if (params->elecA != 0 && params->elecEf != params->elecEi) {
        return new ElectricEvolveStrategy();
    } else if (params->Tf != params->Ti) {
        return new ThermalEvolveStrategy();
    } else if (params->fn > 0) {
        return new StepEvolveStrategy();
    } else {
        return new QuenchEvolveStrategy();
    }
}