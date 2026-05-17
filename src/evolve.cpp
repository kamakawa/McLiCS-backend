#include "../include/evolve.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_eigen.h>
#include <math.h>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "../include/define.h"
#include "../include/geometry.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

// ---------------------------------------------------------------------------
// Internal helper: populate one nni slot from the ni/pt arrays
// ---------------------------------------------------------------------------
static inline void setNni(uint pos, int conditional, nni* nLocal,
                           float* ni, int* pt) {
  if (conditional == 1) {
    nLocal->x  = ni[pos * 3 + 0];
    nLocal->y  = ni[pos * 3 + 1];
    nLocal->z  = ni[pos * 3 + 2];
    nLocal->pt = pt[pos];
  } else {
    nLocal->pt = 0;
  }
}

// ---------------------------------------------------------------------------
// EvolveN — constructor / destructor
// ---------------------------------------------------------------------------
EvolveN::EvolveN(float* ni, int* pt, Parameters* params)
    : ni(ni), pt(pt), params(params),
      Nx(params->Nx), Ny(params->Ny), Nz(params->Nz) {
  // Allocate GSL eigen workspace once; reused by measure_block every MCS step.
  eigen_ws_   = gsl_eigen_symmv_alloc(3);
  eigen_eval_ = gsl_vector_alloc(3);
  eigen_evec_ = gsl_matrix_alloc(3, 3);
  eigen_m_    = gsl_matrix_alloc(3, 3);
}

EvolveN::~EvolveN() {
  gsl_eigen_symmv_free(eigen_ws_);
  gsl_vector_free(eigen_eval_);
  gsl_matrix_free(eigen_evec_);
  gsl_matrix_free(eigen_m_);
}

// ---------------------------------------------------------------------------
// energy_calculator
// ---------------------------------------------------------------------------
float EvolveN::energy_calculator() {
  const int Nt = Nx * Ny * Nz;

  // Count valid (bulk) points once on first call.
  static int  valid      = 0;
  static bool initialized = false;
  if (!initialized) {
    for (int ii = 0; ii < Nt; ii++)
      if (pt[ii]) valid++;
    initialized = true;
  }

  double Etot = 0;
#pragma omp parallel for simd reduction(+:Etot) schedule(simd:dynamic) \
    num_threads(omp_get_max_threads())
  for (int idx = 0; idx < Nt; idx++) {
    nni nLocal[28];

    const int i = idx % Nx;
    const int j = (idx / Nx) % Ny;
    const int k = idx / (Nx * Ny);

    int im = i - 1, ip = i + 1;
    int jm = j - 1, jp = j + 1;
    int km = k - 1, kp = k + 1;
    int neighbour[6] = {
        params->XBound(ip, Nx), params->XBound(im, Nx),
        params->YBound(jp, Ny), params->YBound(jm, Ny),
        params->ZBound(kp, Nz), params->ZBound(km, Nz)};

    setNni(i  + Nx * (j  + Ny * k ), 1,            &nLocal[0], ni, pt);
    setNni(ip + Nx * (j  + Ny * k ), neighbour[0],  &nLocal[1], ni, pt);
    setNni(im + Nx * (j  + Ny * k ), neighbour[1],  &nLocal[2], ni, pt);
    setNni(i  + Nx * (jp + Ny * k ), neighbour[2],  &nLocal[3], ni, pt);
    setNni(i  + Nx * (jm + Ny * k ), neighbour[3],  &nLocal[4], ni, pt);
    setNni(i  + Nx * (j  + Ny * kp), neighbour[4],  &nLocal[5], ni, pt);
    setNni(i  + Nx * (j  + Ny * km), neighbour[5],  &nLocal[6], ni, pt);

    if (nLocal[0].pt > 1) {
      const int base = (i + Nx * (j + Ny * k)) * 3;
      nLocal[7] = { geometry->ns[base + 0],
                    geometry->ns[base + 1],
                    geometry->ns[base + 2], 1 };
    }

    if (params->neighbourKind > 1) {
      setNni(ip + Nx * (jp + Ny * k ), neighbour[0]*neighbour[2], &nLocal[ 8], ni, pt);
      setNni(ip + Nx * (jm + Ny * k ), neighbour[0]*neighbour[3], &nLocal[ 9], ni, pt);
      setNni(ip + Nx * (j  + Ny * kp), neighbour[0]*neighbour[4], &nLocal[10], ni, pt);
      setNni(ip + Nx * (j  + Ny * km), neighbour[0]*neighbour[5], &nLocal[11], ni, pt);
      setNni(im + Nx * (jp + Ny * k ), neighbour[1]*neighbour[2], &nLocal[12], ni, pt);
      setNni(im + Nx * (jm + Ny * k ), neighbour[1]*neighbour[3], &nLocal[13], ni, pt);
      setNni(im + Nx * (j  + Ny * kp), neighbour[1]*neighbour[4], &nLocal[14], ni, pt);
      setNni(im + Nx * (j  + Ny * km), neighbour[1]*neighbour[5], &nLocal[15], ni, pt);
      setNni(i  + Nx * (jp + Ny * kp), neighbour[2]*neighbour[4], &nLocal[16], ni, pt);
      setNni(i  + Nx * (jp + Ny * km), neighbour[2]*neighbour[5], &nLocal[17], ni, pt);
      setNni(i  + Nx * (jm + Ny * kp), neighbour[3]*neighbour[4], &nLocal[18], ni, pt);
      setNni(i  + Nx * (jm + Ny * km), neighbour[3]*neighbour[5], &nLocal[19], ni, pt);
    }
    if (params->neighbourKind == 3) {
      setNni(ip + Nx*(jp + Ny*kp), neighbour[0]*neighbour[2]*neighbour[4], &nLocal[20], ni, pt);
      setNni(ip + Nx*(jp + Ny*km), neighbour[0]*neighbour[2]*neighbour[5], &nLocal[21], ni, pt);
      setNni(ip + Nx*(jm + Ny*kp), neighbour[0]*neighbour[3]*neighbour[4], &nLocal[22], ni, pt);
      setNni(ip + Nx*(jm + Ny*km), neighbour[0]*neighbour[3]*neighbour[5], &nLocal[23], ni, pt);
      setNni(im + Nx*(jp + Ny*kp), neighbour[1]*neighbour[2]*neighbour[4], &nLocal[24], ni, pt);
      setNni(im + Nx*(jp + Ny*km), neighbour[1]*neighbour[2]*neighbour[5], &nLocal[25], ni, pt);
      setNni(im + Nx*(jm + Ny*kp), neighbour[1]*neighbour[3]*neighbour[4], &nLocal[26], ni, pt);
      setNni(im + Nx*(jm + Ny*km), neighbour[1]*neighbour[3]*neighbour[5], &nLocal[27], ni, pt);
    }

    Etot += geometry->latice_Potential(nLocal);
  }
  return static_cast<float>(Etot / Nt);
}

// ---------------------------------------------------------------------------
// Monte_Carlo_Step
// ---------------------------------------------------------------------------
void EvolveN::Monte_Carlo_Step(float& ang_var, gsl_rng** r) {
  const float T  = params->T;
  const int   Nt = Nx * Ny * Nz;

  static int  valid       = 0;
  static bool initialized = false;
  if (!initialized) {
    for (int ii = 0; ii < Nt; ii++)
      if (pt[ii]) valid++;
    initialized = true;
  }

  int acceptance = 0;
  const int iBoxSize = 2, jBoxSize = 2, kBoxSize = 2;
  const int Ni = Nx / 2 + Nx % 2;
  const int Nj = Ny / 2 + Ny % 2;
  const int Nk = Nz / 2 + Nz % 2;
  const int Ntt = Ni * Nj * Nk;

  const int nti = (params->neighbourKind == 2) ? 20 :
                  (params->neighbourKind == 3) ? 28 : 8;

#pragma omp parallel num_threads(omp_get_max_threads())
  for (int tick = 0; tick < 8; tick++) {
#pragma omp for simd reduction(+:acceptance) schedule(simd:dynamic)
    for (int nt = 0; nt < Ntt; nt++) {
      float E_new, E_old, rotation_type, va, ranVal;
      float nNew[3];
      nni nLocal[28];

      const int thread = omp_get_thread_num();
      const int iBox = nt % Ni;
      const int jBox = (nt / Ni) % Nj;
      const int kBox = nt / Ni / Nj;
      const int di = tick % 2;
      const int dj = (tick / 2) % 2;
      const int dk = tick / 4;

      const int i = di + iBox * iBoxSize; if (i >= Nx) continue;
      const int j = dj + jBox * jBoxSize; if (j >= Ny) continue;
      const int k = dk + kBox * kBoxSize; if (k >= Nz) continue;
      if (pti(i, j, k) == 0) continue;

      int im = i-1, ip = i+1;
      int jm = j-1, jp = j+1;
      int km = k-1, kp = k+1;
      int neighbour[6] = {
          params->XBound(ip, Nx), params->XBound(im, Nx),
          params->YBound(jp, Ny), params->YBound(jm, Ny),
          params->ZBound(kp, Nz), params->ZBound(km, Nz)};

      setNni(i  + Nx*(j  + Ny*k ), 1,            &nLocal[0], ni, pt);
      setNni(ip + Nx*(j  + Ny*k ), neighbour[0],  &nLocal[1], ni, pt);
      setNni(im + Nx*(j  + Ny*k ), neighbour[1],  &nLocal[2], ni, pt);
      setNni(i  + Nx*(jp + Ny*k ), neighbour[2],  &nLocal[3], ni, pt);
      setNni(i  + Nx*(jm + Ny*k ), neighbour[3],  &nLocal[4], ni, pt);
      setNni(i  + Nx*(j  + Ny*kp), neighbour[4],  &nLocal[5], ni, pt);
      setNni(i  + Nx*(j  + Ny*km), neighbour[5],  &nLocal[6], ni, pt);

      if (nLocal[0].pt > 1) {
        const int base = (i + Nx*(j + Ny*k)) * 3;
        nLocal[7] = { geometry->ns[base+0], geometry->ns[base+1],
                      geometry->ns[base+2], 1 };
      }

      if (params->neighbourKind > 1) {
        setNni(ip+Nx*(jp+Ny*k ), neighbour[0]*neighbour[2], &nLocal[ 8], ni, pt);
        setNni(ip+Nx*(jm+Ny*k ), neighbour[0]*neighbour[3], &nLocal[ 9], ni, pt);
        setNni(ip+Nx*(j +Ny*kp), neighbour[0]*neighbour[4], &nLocal[10], ni, pt);
        setNni(ip+Nx*(j +Ny*km), neighbour[0]*neighbour[5], &nLocal[11], ni, pt);
        setNni(im+Nx*(jp+Ny*k ), neighbour[1]*neighbour[2], &nLocal[12], ni, pt);
        setNni(im+Nx*(jm+Ny*k ), neighbour[1]*neighbour[3], &nLocal[13], ni, pt);
        setNni(im+Nx*(j +Ny*kp), neighbour[1]*neighbour[4], &nLocal[14], ni, pt);
        setNni(im+Nx*(j +Ny*km), neighbour[1]*neighbour[5], &nLocal[15], ni, pt);
        setNni(i +Nx*(jp+Ny*kp), neighbour[2]*neighbour[4], &nLocal[16], ni, pt);
        setNni(i +Nx*(jp+Ny*km), neighbour[2]*neighbour[5], &nLocal[17], ni, pt);
        setNni(i +Nx*(jm+Ny*kp), neighbour[3]*neighbour[4], &nLocal[18], ni, pt);
        setNni(i +Nx*(jm+Ny*km), neighbour[3]*neighbour[5], &nLocal[19], ni, pt);
      }
      if (params->neighbourKind == 3) {
        setNni(ip+Nx*(jp+Ny*kp), neighbour[0]*neighbour[2]*neighbour[4], &nLocal[20], ni, pt);
        setNni(ip+Nx*(jp+Ny*km), neighbour[0]*neighbour[2]*neighbour[5], &nLocal[21], ni, pt);
        setNni(ip+Nx*(jm+Ny*kp), neighbour[0]*neighbour[3]*neighbour[4], &nLocal[22], ni, pt);
        setNni(ip+Nx*(jm+Ny*km), neighbour[0]*neighbour[3]*neighbour[5], &nLocal[23], ni, pt);
        setNni(im+Nx*(jp+Ny*kp), neighbour[1]*neighbour[2]*neighbour[4], &nLocal[24], ni, pt);
        setNni(im+Nx*(jp+Ny*km), neighbour[1]*neighbour[2]*neighbour[5], &nLocal[25], ni, pt);
        setNni(im+Nx*(jm+Ny*kp), neighbour[1]*neighbour[3]*neighbour[4], &nLocal[26], ni, pt);
        setNni(im+Nx*(jm+Ny*km), neighbour[1]*neighbour[3]*neighbour[5], &nLocal[27], ni, pt);
      }

      va = (2.f * gsl_rng_uniform(r[thread]) - 1.f) * ang_var;
      rotation_type = gsl_rng_uniform(r[thread]);

      if (rotation_type < 0.333f) {
        nNew[0] = nix(i,j,k);
        nNew[1] = niy(i,j,k)*cosf(M_PI*va) + niz(i,j,k)*sinf(M_PI*va);
        nNew[2] = niz(i,j,k)*cosf(M_PI*va) - niy(i,j,k)*sinf(M_PI*va);
      } else if (rotation_type < 0.666f) {
        nNew[0] = nix(i,j,k)*cosf(M_PI*va) - niz(i,j,k)*sinf(M_PI*va);
        nNew[1] = niy(i,j,k);
        nNew[2] = niz(i,j,k)*cosf(M_PI*va) + nix(i,j,k)*sinf(M_PI*va);
      } else {
        nNew[0] = nix(i,j,k)*cosf(M_PI*va) + niy(i,j,k)*sinf(M_PI*va);
        nNew[1] = niy(i,j,k)*cosf(M_PI*va) - nix(i,j,k)*sinf(M_PI*va);
        nNew[2] = niz(i,j,k);
      }

      E_old = geometry->latice_Potential(nLocal);
      nLocal[0].x = nNew[0];
      nLocal[0].y = nNew[1];
      nLocal[0].z = nNew[2];
      E_new = geometry->latice_Potential(nLocal);

      if (gsl_rng_uniform(r[thread]) < expf(-(E_new - E_old) / T)) {
        nix(i,j,k) = nNew[0];
        niy(i,j,k) = nNew[1];
        niz(i,j,k) = nNew[2];
        acceptance++;
      }
    }
  }

  const float ratio = static_cast<float>(acceptance) / valid;
  if (ratio < 0.5f) {
    ang_var *= 0.99f;
    if (ang_var < 0.01f) ang_var = 0.01f;
  } else {
    ang_var /= 0.99f;
    if (ang_var > 1.0f) ang_var -= 0.5f;
  }
}

// ---------------------------------------------------------------------------
// check_Points
// ---------------------------------------------------------------------------
void Evolve::check_Points(int* pt, Parameters params) {
  VallidPoints = 0;
  const int Nx = params.Nx, Ny = params.Ny, Nz = params.Nz;
  for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
      for (int i = 0; i < Nx; i++)
        if (pti(i, j, k)) VallidPoints++;
}

// ---------------------------------------------------------------------------
// Template Method helpers
// ---------------------------------------------------------------------------

// Run MCT equilibration steps (no measurements).
void EvolveN::equilibrate(int mct, float& ang_var, RngPool& rng) {
  for (int step = 0; step < mct; step++)
    Monte_Carlo_Step(ang_var, rng.data());
}

// Run MCS measurement steps and return averaged Measurement.
// Uses the eigen workspace stored in the object — no static locals.
Measurement EvolveN::measure_block(int mcs, float& ang_var, RngPool& rng) {
  Measurement m;
  float mat_n[9], vec_n[3];

  for (int step = 0; step < mcs; step++) {
    Monte_Carlo_Step(ang_var, rng.data());

    const float e = energy_calculator();
    m.E  += e;
    m.E2 += e * e;

    // Build order-parameter matrix Q using per-object workspace
    Matrice_constructor(ni, mat_n, pt, *params);

    // Fill GSL matrix
    for (int ii = 0; ii < 3; ii++)
      for (int jj = 0; jj < 3; jj++)
        gsl_matrix_set(eigen_m_, jj, ii, mat_n[3*ii + jj]);

    gsl_eigen_symmv(eigen_m_, eigen_eval_, eigen_evec_, eigen_ws_);
    gsl_eigen_symmv_sort(eigen_eval_, eigen_evec_, GSL_EIGEN_SORT_VAL_DESC);
    for (int ii = 0; ii < 3; ii++)
      vec_n[ii] = static_cast<float>(gsl_matrix_get(eigen_evec_, ii, 0));

    const float s = static_cast<float>(gsl_vector_get(eigen_eval_, 0));
    m.S1 += s;
    m.S2 += s * s;
  }

  m.E  /= mcs;  m.E2  /= mcs;
  m.S1 /= mcs;  m.S2  /= mcs;
  return m;
}

// Open po.dat in append mode; write header only if the file is empty.
FILE* EvolveN::open_po_file(const char* header) {
  FILE* f = fopen("po.dat", "a");
  if (!f) { perror("po.dat"); exit(1); }
  if (ftell(f) == 0) fprintf(f, "%s\n", header);
  fflush(f);
  return f;
}

// Write one director-field snapshot.
void EvolveN::save_snapshot(const char* fname) {
  print_n(fname, ni, *params, pt);
}

// Append one measurement row to po.dat; key is T or E or file index.
void EvolveN::log_measurement(FILE* po, float key, const Measurement& m) {
  fprintf(po, "%g %g %g %g %g\n",
          key, m.S1, m.S2 - m.S1*m.S1, m.E, m.E2 - m.E*m.E);
  fflush(po);
}
