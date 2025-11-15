#include "../include/evolve.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <omp.h>
#include <string.h>

#include <iostream>
#include <vector>

#include "../include/define.h"
#include "../include/geometry.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

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

float EvolveN::energy_calculator() {
  float *ni_ptr = ni;
  int *pt_ptr = pt;
  float *ns_ptr = geometry->ns.get(); // ns é unique_ptr<float[]>

  static const int Nt = Nx * Ny * Nz;
  static int valid = 0;
  static bool setUp = true;
  if (setUp) {
    for (int ii = 0; ii < Nt; ii++) {
      if (pt_ptr[ii])
        valid++;
    }
    setUp = false;
  }

  double Etot = 0;
#pragma omp parallel for simd reduction(+ : Etot) schedule(simd : dynamic) num_threads(omp_get_max_threads())
  for (int idx = 0; idx < Nt; idx++) {
    nni nLocal[28];

    int i = idx % Nx;
    int j = (idx / Nx) % Ny;
    int k = idx / (Nx * Ny);

    if (i >= Nx || i<0) {
      printf("Tem cachorro nesse mato (%d %d %d )\n", i, j, k);
      fflush(stdout);
      exit(1);
    }
    if (j >= Ny|| j<0) {
      printf("Tem cachorro nesse mato (%d %d %d )\n", i, j, k);
      fflush(stdout);
      exit(1);
    }
    if (k >= Nz|| k<0) {
      printf("Tem cachorro nesse mato (%d %d %d )\n", i, j, k);
      fflush(stdout);
      exit(1);
    }

    int im = (i - 1), ip = (i + 1);
    int jm = (j - 1), jp = (j + 1);
    int km = (k - 1), kp = (k + 1);
    int neighbour[6] = {
        params->XBound(ip, Nx), params->XBound(im, Nx),
        params->YBound(jp, Ny), params->YBound(jm, Ny),
        params->ZBound(kp, Nz), params->ZBound(km, Nz)};
        
    int pos0 = i + Nx * (j + Ny * k );

    setNni(pos0, 1, &nLocal[0], ni_ptr, pt_ptr);
    setNni(ip+ Nx * (j + Ny * k ), neighbour[0], &nLocal[1], ni_ptr, pt_ptr);
    setNni(im+ Nx * (j + Ny * k ), neighbour[1], &nLocal[2], ni_ptr, pt_ptr);
    setNni(i + Nx * (jp+ Ny * k ), neighbour[2], &nLocal[3], ni_ptr, pt_ptr);
    setNni(i + Nx * (jm+ Ny * k ), neighbour[3], &nLocal[4], ni_ptr, pt_ptr);
    setNni(i + Nx * (j + Ny * kp), neighbour[4], &nLocal[5], ni_ptr, pt_ptr);
    setNni(i + Nx * (j + Ny * km), neighbour[5], &nLocal[6], ni_ptr, pt_ptr);

    if (nLocal[0].pt > 1) {
      nLocal[7].x = ns_ptr[pos0 * 3 + 0];
      nLocal[7].y = ns_ptr[pos0 * 3 + 1];
      nLocal[7].z = ns_ptr[pos0 * 3 + 2];
      nLocal[7].pt = 1;
    }
    if (params->neighbourKind > 1) {
      setNni(ip+ Nx * (jp+ Ny * k ), neighbour[0] * neighbour[2], &nLocal[ 8], ni_ptr, pt_ptr);
      setNni(ip+ Nx * (jm+ Ny * k ), neighbour[0] * neighbour[3], &nLocal[ 9], ni_ptr, pt_ptr);
      setNni(ip+ Nx * (j + Ny * kp), neighbour[0] * neighbour[4], &nLocal[10], ni_ptr, pt_ptr);
      setNni(ip+ Nx * (j + Ny * km), neighbour[0] * neighbour[5], &nLocal[11], ni_ptr, pt_ptr);
      setNni(im+ Nx * (jp+ Ny * k ), neighbour[1] * neighbour[2], &nLocal[12], ni_ptr, pt_ptr);
      setNni(im+ Nx * (jm+ Ny * k ), neighbour[1] * neighbour[3], &nLocal[13], ni_ptr, pt_ptr);
      setNni(im+ Nx * (j + Ny * kp), neighbour[1] * neighbour[4], &nLocal[14], ni_ptr, pt_ptr);
      setNni(im+ Nx * (j + Ny * km), neighbour[1] * neighbour[5], &nLocal[15], ni_ptr, pt_ptr);
      setNni(i + Nx * (jp+ Ny * kp), neighbour[2] * neighbour[4], &nLocal[16], ni_ptr, pt_ptr);
      setNni(i + Nx * (jp+ Ny * km), neighbour[2] * neighbour[5], &nLocal[17], ni_ptr, pt_ptr);
      setNni(i + Nx * (jm+ Ny * kp), neighbour[3] * neighbour[4], &nLocal[18], ni_ptr, pt_ptr);
      setNni(i + Nx * (jm+ Ny * km), neighbour[3] * neighbour[5], &nLocal[19], ni_ptr, pt_ptr);
    }
    if (params->neighbourKind == 3) {
      setNni(ip + Nx * (jp + Ny * kp), neighbour[0] * neighbour[2] * neighbour[4], &nLocal[20], ni_ptr, pt_ptr);
      setNni(ip + Nx * (jp + Ny * km), neighbour[0] * neighbour[2] * neighbour[5], &nLocal[21], ni_ptr, pt_ptr);
      setNni(ip + Nx * (jm + Ny * kp), neighbour[0] * neighbour[3] * neighbour[4], &nLocal[22], ni_ptr, pt_ptr);
      setNni(ip + Nx * (jm + Ny * km), neighbour[0] * neighbour[3] * neighbour[5], &nLocal[23], ni_ptr, pt_ptr);
      setNni(im + Nx * (jp + Ny * kp), neighbour[1] * neighbour[2] * neighbour[4], &nLocal[24], ni_ptr, pt_ptr);
      setNni(im + Nx * (jp + Ny * km), neighbour[1] * neighbour[2] * neighbour[5], &nLocal[25], ni_ptr, pt_ptr);
      setNni(im + Nx * (jm + Ny * kp), neighbour[1] * neighbour[3] * neighbour[4], &nLocal[26], ni_ptr, pt_ptr);
      setNni(im + Nx * (jm + Ny * km), neighbour[1] * neighbour[3] * neighbour[5], &nLocal[27], ni_ptr, pt_ptr);
    }
    Etot += geometry->lattice_Potential(nLocal);
  }
  return Etot / Nt;
}

void EvolveN::Monte_Carlo_Step(float &ang_var, gsl_rng **r) {
  float T = params->T;
  static int Nt = Nx * Ny * Nz;
  static int vallid = 0;

  std::vector<float> myNi(Nt * 3);
  
  int acceptance;
  static int time = 0;
  acceptance = 0;
  int iBoxSize = 2;
  int jBoxSize = 2;
  int kBoxSize = 2;
  int Ni = Nx / 2 + Nx % 2;
  int Nj = Ny / 2 + Ny % 2;
  int Nk = Nz / 2 + Nz % 2;
  int Ntt = Ni * Nj * Nk;
  static bool setUp = true;

  float *ni_ptr = ni;
  int *pt_ptr = pt;
  float *ns_ptr = geometry->ns.get(); 

  if (setUp) {
    for (int ii = 0; ii < Nt; ii++)
      if (pt_ptr[ii])
        vallid++;
    setUp = false;
  }
  
  std::copy(ni_ptr, ni_ptr + Nt * 3, myNi.begin());

#pragma omp parallel num_threads(omp_get_max_threads())
  for (int tick = 0; tick < 8; tick++) {
#pragma omp for simd reduction(+ : acceptance) schedule(simd : dynamic)
    for (int nt = 0; nt < Ntt; nt++) {
      float E_new, E_old, rotation_type, va;
      int i, j, k;
      float nNew[3];
      int thread = omp_get_thread_num();
      int iBox = nt % Ni;
      int jBox = (nt / Ni) % Nj;
      int kBox = nt / Ni / Nj;
      int di = tick % 2;
      int dj = (tick / 2) % 2;
      int dk = tick / 4;
      const static int nti = params->neighbourKind == 2 ? 20 : (params->neighbourKind == 3 ? +28 : 8);
      nni nLocal[nti];

      i = di + iBox * iBoxSize;
      if (i >= Nx)
        continue;
      j = dj + jBox * jBoxSize;
      if (j >= Ny)
        continue;
      k = dk + kBox * kBoxSize;
      if (k >= Nz)
        continue;
      if (pt_ptr[i + Nx * (j + Ny * k)] == 0)
        continue;
        
      int im = (i - 1), ip = (i + 1);
      int jm = (j - 1), jp = (j + 1);
      int km = (k - 1), kp = (k + 1);
      int neighbour[6] = {
          params->XBound(ip, Nx), params->XBound(im, Nx),
          params->YBound(jp, Ny), params->YBound(jm, Ny),
          params->ZBound(kp, Nz), params->ZBound(km, Nz)};
          
      int pos0 = i + Nx * (j + Ny * k );
          
      setNni(pos0, 1, &nLocal[0], ni_ptr, pt_ptr);
      setNni(ip+ Nx * (j + Ny * k ), neighbour[0], &nLocal[1], ni_ptr, pt_ptr);
      setNni(im+ Nx * (j + Ny * k ), neighbour[1], &nLocal[2], ni_ptr, pt_ptr);
      setNni(i + Nx * (jp+ Ny * k ), neighbour[2], &nLocal[3], ni_ptr, pt_ptr);
      setNni(i + Nx * (jm+ Ny * k ), neighbour[3], &nLocal[4], ni_ptr, pt_ptr);
      setNni(i + Nx * (j + Ny * kp), neighbour[4], &nLocal[5], ni_ptr, pt_ptr);
      setNni(i + Nx * (j + Ny * km), neighbour[5], &nLocal[6], ni_ptr, pt_ptr);
      
      if (nLocal[0].pt > 1) {
        nLocal[7].x = ns_ptr[pos0 * 3 + 0];
        nLocal[7].y = ns_ptr[pos0 * 3 + 1];
        nLocal[7].z = ns_ptr[pos0 * 3 + 2];
        nLocal[7].pt = 1;
      }
      
      if (params->neighbourKind > 1) {
        setNni(ip+ Nx * (jp+ Ny * k ), neighbour[0] * neighbour[2], &nLocal[ 8], ni_ptr, pt_ptr);
        setNni(ip+ Nx * (jm+ Ny * k ), neighbour[0] * neighbour[3], &nLocal[ 9], ni_ptr, pt_ptr);
        setNni(ip+ Nx * (j + Ny * kp), neighbour[0] * neighbour[4], &nLocal[10], ni_ptr, pt_ptr);
        setNni(ip+ Nx * (j + Ny * km), neighbour[0] * neighbour[5], &nLocal[11], ni_ptr, pt_ptr);
        setNni(im+ Nx * (jp+ Ny * k ), neighbour[1] * neighbour[2], &nLocal[12], ni_ptr, pt_ptr);
        setNni(im+ Nx * (jm+ Ny * k ), neighbour[1] * neighbour[3], &nLocal[13], ni_ptr, pt_ptr);
        setNni(im+ Nx * (j + Ny * kp), neighbour[1] * neighbour[4], &nLocal[14], ni_ptr, pt_ptr);
        setNni(im+ Nx * (j + Ny * km), neighbour[1] * neighbour[5], &nLocal[15], ni_ptr, pt_ptr);
        setNni(i + Nx * (jp+ Ny * kp), neighbour[2] * neighbour[4], &nLocal[16], ni_ptr, pt_ptr);
        setNni(i + Nx * (jp+ Ny * km), neighbour[2] * neighbour[5], &nLocal[17], ni_ptr, pt_ptr);
        setNni(i + Nx * (jm+ Ny * kp), neighbour[3] * neighbour[4], &nLocal[18], ni_ptr, pt_ptr);
        setNni(i + Nx * (jm+ Ny * km), neighbour[3] * neighbour[5], &nLocal[19], ni_ptr, pt_ptr);
      }
      if (params->neighbourKind == 3) {
        setNni(ip + Nx * (jp + Ny * kp), neighbour[0] * neighbour[2] * neighbour[4], &nLocal[20], ni_ptr, pt_ptr);
        setNni(ip + Nx * (jp + Ny * km), neighbour[0] * neighbour[2] * neighbour[5], &nLocal[21], ni_ptr, pt_ptr);
        setNni(ip + Nx * (jm + Ny * kp), neighbour[0] * neighbour[3] * neighbour[4], &nLocal[22], ni_ptr, pt_ptr);
        setNni(ip + Nx * (jm + Ny * km), neighbour[0] * neighbour[3] * neighbour[5], &nLocal[23], ni_ptr, pt_ptr);
        setNni(im + Nx * (jp + Ny * kp), neighbour[1] * neighbour[2] * neighbour[4], &nLocal[24], ni_ptr, pt_ptr);
        setNni(im + Nx * (jp + Ny * km), neighbour[1] * neighbour[2] * neighbour[5], &nLocal[25], ni_ptr, pt_ptr);
        setNni(im + Nx * (jm + Ny * kp), neighbour[1] * neighbour[3] * neighbour[4], &nLocal[26], ni_ptr, pt_ptr);
        setNni(im + Nx * (jm + Ny * km), neighbour[1] * neighbour[3] * neighbour[5], &nLocal[27], ni_ptr, pt_ptr);
      }
      
      va = (2 * gsl_rng_uniform(r[thread]) - 1) * ang_var;
      rotation_type = 0.7;
      {
        int idx_local = pos0 * 3;
        float nx_old = ni_ptr[idx_local + 0];
        float ny_old = ni_ptr[idx_local + 1];
        float nz_old = ni_ptr[idx_local + 2];
        
        if (rotation_type < 0.333) {
          nNew[0] = nx_old;
          nNew[1] = ny_old * cosf(M_PI * va) + nz_old * sinf(M_PI * va);
          nNew[2] = nz_old * cosf(M_PI * va) - ny_old * sinf(M_PI * va);
        } else if (rotation_type < 0.666) {
          nNew[0] = nx_old * cosf(M_PI * va) - nz_old * sinf(M_PI * va);
          nNew[1] = ny_old;
          nNew[2] = nz_old * cosf(M_PI * va) + nx_old * sinf(M_PI * va);
        } else {
          nNew[0] = nx_old * cosf(M_PI * va) + ny_old * sinf(M_PI * va);
          nNew[1] = ny_old * cosf(M_PI * va) - nx_old * sinf(M_PI * va);
          nNew[2] = nz_old;
        }
      }

      E_old = geometry->lattice_Potential(nLocal);
      
      nLocal[0].x = nNew[0];
      nLocal[0].y = nNew[1];
      nLocal[0].z = nNew[2];
      E_new = geometry->lattice_Potential(nLocal);

      if (gsl_rng_uniform(r[thread]) < exp(-(E_new - E_old) / T)) {
        int idx_local = pos0 * 3;
        ni_ptr[idx_local + 0] = nNew[0];
        ni_ptr[idx_local + 1] = nNew[1];
        ni_ptr[idx_local + 2] = nNew[2];
        acceptance++;
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
      ang_var -= 0.5;
    }
  }
}