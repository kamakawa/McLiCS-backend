#include "../../include/cuda/mc_bulk_step.cuh"

#ifdef USE_CUDA

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdint>

#include "../../include/cuda/potential_device.cuh"
#include "../../include/cuda/params_device.cuh"

#ifndef CUDA_CHECK
#define CUDA_CHECK(call) do { \
  cudaError_t _e = (call); \
  if (_e != cudaSuccess) { \
    std::fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(_e)); \
    std::abort(); \
  } \
} while(0)
#endif

// ---------------- RNG simples por thread (xorshift32) ----------------
__device__ __forceinline__ uint32_t xorshift32(uint32_t& s) {
  uint32_t x = s;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  s = x;
  return x;
}
__device__ __forceinline__ float u01(uint32_t& s) {
  // 24 bits
  return (xorshift32(s) >> 8) * (1.0f / 16777216.0f);
}

// ---------------- splitmix64 para espalhar seed ----------------
__device__ __forceinline__ std::uint64_t splitmix64_device(std::uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}
__device__ __forceinline__ uint32_t splitmix32_device(std::uint64_t x) {
  return (uint32_t)(splitmix64_device(x) >> 32);
}

// ---------------- rotação (igual ao CPU) ----------------
__device__ __forceinline__ void rotate_vector(float x, float y, float z,
                                              float va, float rotType,
                                              float& ox, float& oy, float& oz) {
  const float c = cosf((float)M_PI * va);
  const float s = sinf((float)M_PI * va);

  if (rotType < 0.33333334f) {
    ox = x;
    oy = y * c + z * s;
    oz = z * c - y * s;
  } else if (rotType < 0.66666667f) {
    ox = x * c - z * s;
    oy = y;
    oz = z * c + x * s;
  } else {
    ox = x * c + y * s;
    oy = y * c - x * s;
    oz = z;
  }

  const float invn = rsqrtf(ox*ox + oy*oy + oz*oz);
  ox *= invn; oy *= invn; oz *= invn;
}

// ---------------- seletor de energia bulk (arquitetura robusta) ----------------
__device__ __forceinline__ float bulk_energy_dispatch(
    BulkPotentialType t,
    const float ni[3], const float nj[3],
    const ParamsDevice* p,
    const float rij[3],
    int nk
) {
  switch (t) {
    case BulkPotentialType::LL:
      return bulk_energy_LL(ni, nj, p);
    case BulkPotentialType::GHRL:
      return bulk_energy_GHRL(ni, nj, p, rij, nk);
    case BulkPotentialType::PEAR:
      return bulk_energy_PEAR(ni, nj, p, rij);
    default:
      return bulk_energy_LL(ni, nj, p);
  }
}

// ---------------- kernel: checkerboard 2x2x2 ----------------
__global__ void mc_bulk_kernel_tick(
    float* ni, const int* pt,
    const ParamsDevice* p,
    float ang_var,
    std::uint64_t seed,
    int tick,
    int* acceptance_out
) {
  const int Nx = p->Nx, Ny = p->Ny, Nz = p->Nz;

  const int Ni = Nx/2 + (Nx%2);
  const int Nj = Ny/2 + (Ny%2);
  const int Nk = Nz/2 + (Nz%2);
  const int Ntt = Ni*Nj*Nk;

  const int nt = (int)(blockIdx.x * blockDim.x + threadIdx.x);
  if (nt >= Ntt) return;

  const int iBox = nt % Ni;
  const int jBox = (nt / Ni) % Nj;
  const int kBox = nt / (Ni * Nj);

  const int di = tick % 2;
  const int dj = (tick / 2) % 2;
  const int dk = tick / 4;

  const int i = di + iBox * 2;
  const int j = dj + jBox * 2;
  const int k = dk + kBox * 2;

  if (i >= Nx || j >= Ny || k >= Nz) return;

  const int idx = i + Nx * (j + Ny * k);
  if (pt[idx] == 0) return;

  // RNG (por thread+tick; seed vem do host)
  // OBS: splitmix espalha os bits do seed e reduz correlação => melhora varE/varS
  const std::uint64_t key =
      seed ^ (0xD1B54A32D192ED03ULL
              + (std::uint64_t)nt   * 0x9E3779B97F4A7C15ULL
              + (std::uint64_t)tick * 0xBF58476D1CE4E5B9ULL);
  uint32_t rng = splitmix32_device(key);

  const float va = (2.0f * u01(rng) - 1.0f) * ang_var;
  const float rotType = u01(rng);

  const int base = idx * 3;
  const float ox = ni[base + 0];
  const float oy = ni[base + 1];
  const float oz = ni[base + 2];

  float nx, ny, nz;
  rotate_vector(ox, oy, oz, va, rotType, nx, ny, nz);

  const float n_old[3] = {ox, oy, oz};
  const float n_new[3] = {nx, ny, nz};

  // >>> CRÍTICO p/ variância: acumular em double <<<
  double Eold = 0.0, Enew = 0.0;

  // ---------------- 1o vizinhos (6) ----------------
  {
    const int dii[6] = {+1,-1, 0, 0, 0, 0};
    const int djj[6] = { 0, 0,+1,-1, 0, 0};
    const int dkk[6] = { 0, 0, 0, 0,+1,-1};

    const float rijv[6][3] = {
      {+1,0,0},{-1,0,0},{0,+1,0},{0,-1,0},{0,0,+1},{0,0,-1}
    };

    for (int t = 0; t < 6; ++t) {
      int ii = i + dii[t];
      int jj = j + djj[t];
      int kk = k + dkk[t];

      int okx = boundary_apply(p->bx, ii, Nx);
      int oky = boundary_apply(p->by, jj, Ny);
      int okz = boundary_apply(p->bz, kk, Nz);
      if (!(okx && oky && okz)) continue;

      const int nidx = ii + Nx * (jj + Ny * kk);
      if (pt[nidx] == 0) continue;

      const int nbase = nidx * 3;
      const float n_nb[3] = {ni[nbase+0], ni[nbase+1], ni[nbase+2]};

      Eold += (double)bulk_energy_dispatch(p->bulkType, n_old, n_nb, p, rijv[t], 1);
      Enew += (double)bulk_energy_dispatch(p->bulkType, n_new, n_nb, p, rijv[t], 1);
    }
  }

  // ---------------- 2o vizinhos (12) ----------------
  if (p->neighbourKind > 1) {
    const float is2 = 0.707106781f;

    const int off[12][3] = {
      {+1,+1,0},{+1,-1,0},{-1,+1,0},{-1,-1,0},
      {+1,0,+1},{+1,0,-1},{-1,0,+1},{-1,0,-1},
      {0,+1,+1},{0,+1,-1},{0,-1,+1},{0,-1,-1}
    };

    const float rij[12][3] = {
      {+is2,+is2,0},{+is2,-is2,0},{-is2,+is2,0},{-is2,-is2,0},
      {+is2,0,+is2},{+is2,0,-is2},{-is2,0,+is2},{-is2,0,-is2},
      {0,+is2,+is2},{0,+is2,-is2},{0,-is2,+is2},{0,-is2,-is2}
    };

    for (int t=0;t<12;++t){
      int ii=i+off[t][0], jj=j+off[t][1], kk=k+off[t][2];
      int okx = boundary_apply(p->bx, ii, Nx);
      int oky = boundary_apply(p->by, jj, Ny);
      int okz = boundary_apply(p->bz, kk, Nz);
      if(!(okx&&oky&&okz)) continue;

      const int nidx = ii + Nx * (jj + Ny * kk);
      if (pt[nidx]==0) continue;

      const int nbase=nidx*3;
      const float n_nb[3] = {ni[nbase+0], ni[nbase+1], ni[nbase+2]};

      Eold += (double)bulk_energy_dispatch(p->bulkType, n_old, n_nb, p, rij[t], 2);
      Enew += (double)bulk_energy_dispatch(p->bulkType, n_new, n_nb, p, rij[t], 2);
    }
  }

  // ---------------- 3o vizinhos (8) ----------------
  if (p->neighbourKind == 3) {
    const float is3 = 0.577350269f;

    const int off[8][3] = {
      {+1,+1,+1},{+1,+1,-1},{+1,-1,+1},{+1,-1,-1},
      {-1,+1,+1},{-1,+1,-1},{-1,-1,+1},{-1,-1,-1}
    };

    const float rij[8][3] = {
      {+is3,+is3,+is3},{+is3,+is3,-is3},{+is3,-is3,+is3},{+is3,-is3,-is3},
      {-is3,+is3,+is3},{-is3,+is3,-is3},{-is3,-is3,+is3},{-is3,-is3,-is3}
    };

    for(int t=0;t<8;++t){
      int ii=i+off[t][0], jj=j+off[t][1], kk=k+off[t][2];
      int okx = boundary_apply(p->bx, ii, Nx);
      int oky = boundary_apply(p->by, jj, Ny);
      int okz = boundary_apply(p->bz, kk, Nz);
      if(!(okx&&oky&&okz)) continue;

      const int nidx = ii + Nx * (jj + Ny * kk);
      if (pt[nidx]==0) continue;

      const int nbase=nidx*3;
      const float n_nb[3] = {ni[nbase+0], ni[nbase+1], ni[nbase+2]};

      Eold += (double)bulk_energy_dispatch(p->bulkType, n_old, n_nb, p, rij[t], 3);
      Enew += (double)bulk_energy_dispatch(p->bulkType, n_new, n_nb, p, rij[t], 3);
    }
  }

  // ---------------- campo elétrico ----------------
  Eold += (double)electric_energy(n_old, p);
  Enew += (double)electric_energy(n_new, p);

  const double dE = Enew - Eold;
  const float u = u01(rng);

  // Metropolis (double no exp => melhora estabilidade estatística)
  if (dE < 0.0 || u < (float)exp(-dE / (double)p->T)) {
    ni[base+0]=nx; ni[base+1]=ny; ni[base+2]=nz;
    atomicAdd(acceptance_out, 1);
  }
}

// ---------------- estado device persistente ----------------
struct DeviceState {
  float* d_ni = nullptr;
  int*   d_pt = nullptr;
  int*   d_acc = nullptr;
  ParamsDevice* d_params = nullptr;

  int Nx=0, Ny=0, Nz=0;
  std::size_t bytes_ni = 0;
  std::size_t bytes_pt = 0;
  bool inited = false;
};

static DeviceState g_state;

static void ensure_buffers(int Nx, int Ny, int Nz) {
  const std::size_t N = (std::size_t)Nx * (std::size_t)Ny * (std::size_t)Nz;

  const std::size_t bytes_ni = 3ULL * N * sizeof(float);
  const std::size_t bytes_pt = 1ULL * N * sizeof(int); // só N (pt[idx]) é usado

  if (g_state.inited && g_state.Nx==Nx && g_state.Ny==Ny && g_state.Nz==Nz) {
    g_state.bytes_ni = bytes_ni;
    g_state.bytes_pt = bytes_pt;
    return;
  }

  if (g_state.d_ni) CUDA_CHECK(cudaFree(g_state.d_ni));
  if (g_state.d_pt) CUDA_CHECK(cudaFree(g_state.d_pt));
  if (g_state.d_acc) CUDA_CHECK(cudaFree(g_state.d_acc));
  if (g_state.d_params) CUDA_CHECK(cudaFree(g_state.d_params));

  g_state.Nx=Nx; g_state.Ny=Ny; g_state.Nz=Nz;
  g_state.bytes_ni = bytes_ni;
  g_state.bytes_pt = bytes_pt;

  CUDA_CHECK(cudaMalloc((void**)&g_state.d_ni, bytes_ni));
  CUDA_CHECK(cudaMalloc((void**)&g_state.d_pt, bytes_pt));
  CUDA_CHECK(cudaMalloc((void**)&g_state.d_acc, sizeof(int)));
  CUDA_CHECK(cudaMalloc((void**)&g_state.d_params, sizeof(ParamsDevice)));

  g_state.inited = true;
}

extern "C" void mc_bulk_step_inplace(
    float* h_ni, int* h_pt,
    const ParamsDevice& hp,
    float ang_var,
    std::uint64_t seed,
    int& acceptance_out
) {
  ensure_buffers(hp.Nx, hp.Ny, hp.Nz);

  // upload params
  CUDA_CHECK(cudaMemcpy(g_state.d_params, &hp, sizeof(ParamsDevice), cudaMemcpyHostToDevice));

  // upload lattice
  CUDA_CHECK(cudaMemcpy(g_state.d_ni, h_ni, g_state.bytes_ni, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(g_state.d_pt, h_pt, g_state.bytes_pt, cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaMemset(g_state.d_acc, 0, sizeof(int)));

  const int Nx=hp.Nx, Ny=hp.Ny, Nz=hp.Nz;
  const int Ni = Nx/2 + (Nx%2);
  const int Nj = Ny/2 + (Ny%2);
  const int Nk = Nz/2 + (Nz%2);
  const int Ntt = Ni*Nj*Nk;

  const int block = 256;
  const int grid  = (Ntt + block - 1) / block;

  for (int tick=0; tick<8; ++tick) {
    mc_bulk_kernel_tick<<<grid, block>>>(
      g_state.d_ni, g_state.d_pt,
      g_state.d_params,
      ang_var, seed, tick,
      g_state.d_acc
    );
    CUDA_CHECK(cudaGetLastError());
  }

  int h_acc=0;
  CUDA_CHECK(cudaMemcpy(&h_acc, g_state.d_acc, sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(h_ni, g_state.d_ni, g_state.bytes_ni, cudaMemcpyDeviceToHost));

  acceptance_out = h_acc;
}

#endif // USE_CUDA