#include "../include/evolve.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <cmath>
#include <omp.h>

#include <iostream>
#include <memory>

#include "../include/define.h"
#include "../include/geometry.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

void inline setNni(uint pos, int conditional, nni *nLocal, std::unique_ptr<float[]>& ni, std::unique_ptr<int[]>& pt, int Nx, int Ny) {
  if (conditional == 1) {
    nLocal->x = ni[pos * 3 + 0];
    nLocal->y = ni[pos * 3 + 1];
    nLocal->z = ni[pos * 3 + 2];
    nLocal->pt = pt[pos];
  } else {
    nLocal->pt = 0;
  }
}

float EvolveN::energy_calculator() {
  const int Nt = Nx * Ny * Nz;
  static int valid = 0;
  static bool setUp = false;
  
  if (!setUp) {
    for (int ii = 0; ii < Nt; ii++) {
      if (pt[ii])
        valid++;
    }
    setUp = true;
  }

  double Etot = 0;
#pragma omp parallel for reduction(+ : Etot) schedule(dynamic) num_threads(omp_get_max_threads())
  for (int idx = 0; idx < Nt; idx++) {
    nni nLocal[28];

    int i = idx % Nx;
    int j = (idx / Nx) % Ny;
    int k = idx / (Nx * Ny);
    
    // Verificação de limites removida - já é feita pelos bounds

    int im = (i - 1), ip = (i + 1);
    int jm = (j - 1), jp = (j + 1);
    int km = (k - 1), kp = (k + 1);
    
    int neighbour[6] = {
        params->XBound(ip, Nx), params->XBound(im, Nx),
        params->YBound(jp, Ny), params->YBound(jm, Ny),
        params->ZBound(kp, Nz), params->ZBound(km, Nz)};
        
    setNni(i + Nx * (j + Ny * k ), 1, &nLocal[0], ni, pt, Nx, Ny);
    setNni(ip + Nx * (j + Ny * k ), neighbour[0], &nLocal[1], ni, pt, Nx, Ny);
    setNni(im + Nx * (j + Ny * k ), neighbour[1], &nLocal[2], ni, pt, Nx, Ny);
    setNni(i + Nx * (jp + Ny * k ), neighbour[2], &nLocal[3], ni, pt, Nx, Ny);
    setNni(i + Nx * (jm + Ny * k ), neighbour[3], &nLocal[4], ni, pt, Nx, Ny);
    setNni(i + Nx * (j + Ny * kp), neighbour[4], &nLocal[5], ni, pt, Nx, Ny);
    setNni(i + Nx * (j + Ny * km), neighbour[5], &nLocal[6], ni, pt, Nx, Ny);
    
    if (nLocal[0].pt > 1) {
      nLocal[7].x = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 0];
      nLocal[7].y = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 1];
      nLocal[7].z = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 2];
      nLocal[7].pt = 1;
    }
    
    if (params->neighbourKind > 1) {
      setNni(ip + Nx * (jp + Ny * k ), neighbour[0] * neighbour[2], &nLocal[ 8], ni, pt, Nx, Ny);
      setNni(ip + Nx * (jm + Ny * k ), neighbour[0] * neighbour[3], &nLocal[ 9], ni, pt, Nx, Ny);
      setNni(ip + Nx * (j + Ny * kp), neighbour[0] * neighbour[4], &nLocal[10], ni, pt, Nx, Ny);
      setNni(ip + Nx * (j + Ny * km), neighbour[0] * neighbour[5], &nLocal[11], ni, pt, Nx, Ny);
      setNni(im + Nx * (jp + Ny * k ), neighbour[1] * neighbour[2], &nLocal[12], ni, pt, Nx, Ny);
      setNni(im + Nx * (jm + Ny * k ), neighbour[1] * neighbour[3], &nLocal[13], ni, pt, Nx, Ny);
      setNni(im + Nx * (j + Ny * kp), neighbour[1] * neighbour[4], &nLocal[14], ni, pt, Nx, Ny);
      setNni(im + Nx * (j + Ny * km), neighbour[1] * neighbour[5], &nLocal[15], ni, pt, Nx, Ny);
      setNni(i + Nx * (jp + Ny * kp), neighbour[2] * neighbour[4], &nLocal[16], ni, pt, Nx, Ny);
      setNni(i + Nx * (jp + Ny * km), neighbour[2] * neighbour[5], &nLocal[17], ni, pt, Nx, Ny);
      setNni(i + Nx * (jm + Ny * kp), neighbour[3] * neighbour[4], &nLocal[18], ni, pt, Nx, Ny);
      setNni(i + Nx * (jm + Ny * km), neighbour[3] * neighbour[5], &nLocal[19], ni, pt, Nx, Ny);
    }
    
    if (params->neighbourKind == 3) {
      setNni(ip + Nx * (jp + Ny * kp), neighbour[0] * neighbour[2] * neighbour[4], &nLocal[20], ni, pt, Nx, Ny);
      setNni(ip + Nx * (jp + Ny * km), neighbour[0] * neighbour[2] * neighbour[5], &nLocal[21], ni, pt, Nx, Ny);
      setNni(ip + Nx * (jm + Ny * kp), neighbour[0] * neighbour[3] * neighbour[4], &nLocal[22], ni, pt, Nx, Ny);
      setNni(ip + Nx * (jm + Ny * km), neighbour[0] * neighbour[3] * neighbour[5], &nLocal[23], ni, pt, Nx, Ny);
      setNni(im + Nx * (jp + Ny * kp), neighbour[1] * neighbour[2] * neighbour[4], &nLocal[24], ni, pt, Nx, Ny);
      setNni(im + Nx * (jp + Ny * km), neighbour[1] * neighbour[2] * neighbour[5], &nLocal[25], ni, pt, Nx, Ny);
      setNni(im + Nx * (jm + Ny * kp), neighbour[1] * neighbour[3] * neighbour[4], &nLocal[26], ni, pt, Nx, Ny);
      setNni(im + Nx * (jm + Ny * km), neighbour[1] * neighbour[3] * neighbour[5], &nLocal[27], ni, pt, Nx, Ny);
    }
    
    Etot += geometry->lattice_Potential(nLocal);
  }
  return Etot / Nt;
}

void EvolveN::Monte_Carlo_Step(float &ang_var, gsl_rng **r) {
  float T = params->T;
  const int Nt = Nx * Ny * Nz;
  static int vallid = 0;
  static bool setUp = true;
  
  if (setUp) {
    for (int ii = 0; ii < Nt; ii++)
      if (pt[ii])
        vallid++;
    setUp = false;
  }

  int acceptance = 0;
  const int iBoxSize = 2;
  const int jBoxSize = 2;
  const int kBoxSize = 2;
  const int Ni = Nx / 2 + Nx % 2;
  const int Nj = Ny / 2 + Ny % 2;
  const int Nk = Nz / 2 + Nz % 2;
  const int Ntt = Ni * Nj * Nk;

#pragma omp parallel num_threads(omp_get_max_threads())
  {
    for (int tick = 0; tick < 8; tick++) {
#pragma omp for reduction(+ : acceptance) schedule(dynamic)
      for (int nt = 0; nt < Ntt; nt++) {
        float E_new, E_old, rotation_type, va;
        float nNew[3];
        int thread = omp_get_thread_num();
        
        int iBox = nt % Ni;
        int jBox = (nt / Ni) % Nj;
        int kBox = nt / (Ni * Nj);
        
        int di = tick % 2;
        int dj = (tick / 2) % 2;
        int dk = tick / 4;
        
        const int nti = (params->neighbourKind == 2) ? 20 : ((params->neighbourKind == 3) ? 28 : 8);
        nni nLocal[nti];

        int i = di + iBox * iBoxSize;
        if (i >= Nx) continue;
        int j = dj + jBox * jBoxSize;
        if (j >= Ny) continue;
        int k = dk + kBox * kBoxSize;
        if (k >= Nz) continue;
        
        if (pt[i + Nx * (j + Ny * k)] == 0)
          continue;

        int im = (i - 1), ip = (i + 1);
        int jm = (j - 1), jp = (j + 1);
        int km = (k - 1), kp = (k + 1);
        
        int neighbour[6] = {
            params->XBound(ip, Nx), params->XBound(im, Nx),
            params->YBound(jp, Ny), params->YBound(jm, Ny),
            params->ZBound(kp, Nz), params->ZBound(km, Nz)};
            
        setNni(i + Nx * (j + Ny * k ), 1, &nLocal[0], ni, pt, Nx, Ny);
        setNni(ip + Nx * (j + Ny * k ), neighbour[0], &nLocal[1], ni, pt, Nx, Ny);
        setNni(im + Nx * (j + Ny * k ), neighbour[1], &nLocal[2], ni, pt, Nx, Ny);
        setNni(i + Nx * (jp + Ny * k ), neighbour[2], &nLocal[3], ni, pt, Nx, Ny);
        setNni(i + Nx * (jm + Ny * k ), neighbour[3], &nLocal[4], ni, pt, Nx, Ny);
        setNni(i + Nx * (j + Ny * kp), neighbour[4], &nLocal[5], ni, pt, Nx, Ny);
        setNni(i + Nx * (j + Ny * km), neighbour[5], &nLocal[6], ni, pt, Nx, Ny);
        
        if (nLocal[0].pt > 1) {
          nLocal[7].x = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 0];
          nLocal[7].y = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 1];
          nLocal[7].z = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 2];
          nLocal[7].pt = 1;
        }
        
        if (params->neighbourKind > 1) {
          // ... (mesmo padrão do energy_calculator)
          setNni(ip + Nx * (jp + Ny * k ), neighbour[0] * neighbour[2], &nLocal[8], ni, pt, Nx, Ny);
          setNni(ip + Nx * (jm + Ny * k ), neighbour[0] * neighbour[3], &nLocal[9], ni, pt, Nx, Ny);
          setNni(ip + Nx * (j + Ny * kp), neighbour[0] * neighbour[4], &nLocal[10], ni, pt, Nx, Ny);
          setNni(ip + Nx * (j + Ny * km), neighbour[0] * neighbour[5], &nLocal[11], ni, pt, Nx, Ny);
          setNni(im + Nx * (jp + Ny * k ), neighbour[1] * neighbour[2], &nLocal[12], ni, pt, Nx, Ny);
          setNni(im + Nx * (jm + Ny * k ), neighbour[1] * neighbour[3], &nLocal[13], ni, pt, Nx, Ny);
          setNni(im + Nx * (j + Ny * kp), neighbour[1] * neighbour[4], &nLocal[14], ni, pt, Nx, Ny);
          setNni(im + Nx * (j + Ny * km), neighbour[1] * neighbour[5], &nLocal[15], ni, pt, Nx, Ny);
          setNni(i + Nx * (jp + Ny * kp), neighbour[2] * neighbour[4], &nLocal[16], ni, pt, Nx, Ny);
          setNni(i + Nx * (jp + Ny * km), neighbour[2] * neighbour[5], &nLocal[17], ni, pt, Nx, Ny);
          setNni(i + Nx * (jm + Ny * kp), neighbour[3] * neighbour[4], &nLocal[18], ni, pt, Nx, Ny);
          setNni(i + Nx * (jm + Ny * km), neighbour[3] * neighbour[5], &nLocal[19], ni, pt, Nx, Ny);
        }
        
        if (params->neighbourKind == 3) {
          // ... (mesmo padrão do energy_calculator)
          setNni(ip + Nx * (jp + Ny * kp), neighbour[0] * neighbour[2] * neighbour[4], &nLocal[20], ni, pt, Nx, Ny);
          setNni(ip + Nx * (jp + Ny * km), neighbour[0] * neighbour[2] * neighbour[5], &nLocal[21], ni, pt, Nx, Ny);
          setNni(ip + Nx * (jm + Ny * kp), neighbour[0] * neighbour[3] * neighbour[4], &nLocal[22], ni, pt, Nx, Ny);
          setNni(ip + Nx * (jm + Ny * km), neighbour[0] * neighbour[3] * neighbour[5], &nLocal[23], ni, pt, Nx, Ny);
          setNni(im + Nx * (jp + Ny * kp), neighbour[1] * neighbour[2] * neighbour[4], &nLocal[24], ni, pt, Nx, Ny);
          setNni(im + Nx * (jp + Ny * km), neighbour[1] * neighbour[2] * neighbour[5], &nLocal[25], ni, pt, Nx, Ny);
          setNni(im + Nx * (jm + Ny * kp), neighbour[1] * neighbour[3] * neighbour[4], &nLocal[26], ni, pt, Nx, Ny);
          setNni(im + Nx * (jm + Ny * km), neighbour[1] * neighbour[3] * neighbour[5], &nLocal[27], ni, pt, Nx, Ny);
        }

        va = (2 * gsl_rng_uniform(r[thread]) - 1) * ang_var;
        rotation_type = 0.7; // Valor fixo para debug

        if (rotation_type < 0.333) {
          nNew[0] = ni[i + Nx * (j + Ny * k) * 3 + 0];
          nNew[1] = ni[i + Nx * (j + Ny * k) * 3 + 1] * std::cos(M_PI * va) + ni[i + Nx * (j + Ny * k) * 3 + 2] * std::sin(M_PI * va);
          nNew[2] = ni[i + Nx * (j + Ny * k) * 3 + 2] * std::cos(M_PI * va) - ni[i + Nx * (j + Ny * k) * 3 + 1] * std::sin(M_PI * va);
        } else if (rotation_type < 0.666) {
          nNew[0] = ni[i + Nx * (j + Ny * k) * 3 + 0] * std::cos(M_PI * va) - ni[i + Nx * (j + Ny * k) * 3 + 2] * std::sin(M_PI * va);
          nNew[1] = ni[i + Nx * (j + Ny * k) * 3 + 1];
          nNew[2] = ni[i + Nx * (j + Ny * k) * 3 + 2] * std::cos(M_PI * va) + ni[i + Nx * (j + Ny * k) * 3 + 0] * std::sin(M_PI * va);
        } else {
          nNew[0] = ni[i + Nx * (j + Ny * k) * 3 + 0] * std::cos(M_PI * va) + ni[i + Nx * (j + Ny * k) * 3 + 1] * std::sin(M_PI * va);
          nNew[1] = ni[i + Nx * (j + Ny * k) * 3 + 1] * std::cos(M_PI * va) - ni[i + Nx * (j + Ny * k) * 3 + 0] * std::sin(M_PI * va);
          nNew[2] = ni[i + Nx * (j + Ny * k) * 3 + 2];
        }

        E_old = geometry->lattice_Potential(nLocal);
        
        nLocal[0].x = nNew[0];
        nLocal[0].y = nNew[1];
        nLocal[0].z = nNew[2];
        E_new = geometry->lattice_Potential(nLocal);

        if (gsl_rng_uniform(r[thread]) < std::exp(-(E_new - E_old) / T)) {
          ni[i + Nx * (j + Ny * k) * 3 + 0] = nNew[0];
          ni[i + Nx * (j + Ny * k) * 3 + 1] = nNew[1];
          ni[i + Nx * (j + Ny * k) * 3 + 2] = nNew[2];
          acceptance++;
        }
      }
    }
  }

  if (1.0 * acceptance / vallid < 0.5) {
    ang_var *= 0.99;
    if (ang_var < 0.01)
      ang_var = 0.01;
  } else if (1.0 * acceptance / vallid > 0.5) {
    ang_var /= 0.99;
    if (ang_var > 1.0) {
      ang_var = 1.0;
    }
  }
}

void Evolve::check_Points(int *pt, Parameters params) {
  VallidPoints = 0;
  int Nx = params.Nx;
  int Ny = params.Ny;
  int Nz = params.Nz;
  for (int k = 0; k < Nz; k++) {
    for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
        if (pt[i + Nx * (j + Ny * k)])
          VallidPoints++;
      }
    }
  }
}