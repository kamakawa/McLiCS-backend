#include "../include/evolve.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <omp.h>

#include <iostream>

#include "../include/define.h"
#include "../include/geometry.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

// Note: A função setNni é uma função auxiliar e não acessa 'params', logo não precisa de refatoração.
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

// Alteração 1: energy_calculator() não precisa de const& na assinatura pois não recebe 'params' diretamente.
float EvolveN::energy_calculator() {
  // Alteração 2: Acesso a dimensões da grade aninhadas em 'lattice'
  static const int Nt = params->lattice.Nx * params->lattice.Ny * params->lattice.Nz;
  static int valid = 0;
  static bool setUp;
  if (setUp) {
    for (int ii = 0; ii < Nt; ii++) {
      if (pt[ii])
        valid++;
    }
  }

  double Etot = 0;
#pragma omp parallel for simd reduction(+ : Etot) schedule(simd : dynamic) num_threads(omp_get_max_threads())
  for (int idx = 0; idx < Nt; idx++) {
    nni nLocal[28];

    // Alteração 3, 4, 5: Acesso a dimensões da grade aninhadas em 'lattice'
    int i = idx % params->lattice.Nx;
    int j = (idx / params->lattice.Nx) % params->lattice.Ny;
    int k = idx / (params->lattice.Nx * params->lattice.Ny);
    
    // Alteração 6, 7, 8: Acesso a dimensões da grade aninhadas em 'lattice'
    if (i >= params->lattice.Nx || i<0) {
      printf("Tem cachorro nesse mato (%d %d %d )\n", i, j, k);
      fflush(stdout);
      exit(1);
    }
    // Alteração 9, 10, 11: Acesso a dimensões da grade aninhadas em 'lattice'
    if (j >= params->lattice.Ny|| j<0) {
      printf("Tem cachorro nesse mato (%d %d %d )\n", i, j, k);
      fflush(stdout);
      exit(1);
    }
    // Alteração 12, 13, 14: Acesso a dimensões da grade aninhadas em 'lattice'
    if (k >= params->lattice.Nz|| k<0) {
      printf("Tem cachorro nesse mato (%d %d %d )\n", i, j, k);
      fflush(stdout);
      exit(1);
    }

    int im = (i - 1), ip = (i + 1);
    int jm = (j - 1), jp = (j + 1);
    int km = (k - 1), kp = (k + 1);
    
    // Alteração 15-20: Acesso a funções de fronteira e dimensões aninhadas em 'lattice'
    int neighbour[6] = {
        params->lattice.XBound(ip, params->lattice.Nx), params->lattice.XBound(im, params->lattice.Nx),
        params->lattice.YBound(jp, params->lattice.Ny), params->lattice.YBound(jm, params->lattice.Ny),
        params->lattice.ZBound(kp, params->lattice.Nz), params->lattice.ZBound(km, params->lattice.Nz)};
        
    setNni(i + Nx * (j + Ny * k ), 1, &nLocal[0], ni, pt);
    setNni(ip+ Nx * (j + Ny * k ), neighbour[0], &nLocal[1], ni, pt);
    setNni(im+ Nx * (j + Ny * k ), neighbour[1], &nLocal[2], ni, pt);
    setNni(i + Nx * (jp+ Ny * k ), neighbour[2], &nLocal[3], ni, pt);
    setNni(i + Nx * (jm+ Ny * k ), neighbour[3], &nLocal[4], ni, pt);
    setNni(i + Nx * (j + Ny * kp), neighbour[4], &nLocal[5], ni, pt);
    setNni(i + Nx * (j + Ny * km), neighbour[5], &nLocal[6], ni, pt);
    if (nLocal[0].pt > 1) {
      nLocal[7].x = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 0];
      nLocal[7].y = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 1];
      nLocal[7].z = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 2];
      nLocal[7].pt = 1;
      //~ printf("%g %g %g %d\n",nLocal[7].x,nLocal[7].y,nLocal[7].z,nLocal[7].pt );
    }
    
    // Alteração 21: Acesso a neighbourKind aninhado em 'neighbourhood'
    if (params->neighbourhood.neighbourKind > 1) {
      setNni(ip+ Nx * (jp+ Ny * k ), neighbour[0] * neighbour[2], &nLocal[ 8], ni, pt);
      setNni(ip+ Nx * (jm+ Ny * k ), neighbour[0] * neighbour[3], &nLocal[ 9], ni, pt);
      setNni(ip+ Nx * (j + Ny * kp), neighbour[0] * neighbour[4], &nLocal[10], ni, pt);
      setNni(ip+ Nx * (j + Ny * km), neighbour[0] * neighbour[5], &nLocal[11], ni, pt);
      setNni(im+ Nx * (jp+ Ny * k ), neighbour[1] * neighbour[2], &nLocal[12], ni, pt);
      setNni(im+ Nx * (jm+ Ny * k ), neighbour[1] * neighbour[3], &nLocal[13], ni, pt);
      setNni(im+ Nx * (j + Ny * kp), neighbour[1] * neighbour[4], &nLocal[14], ni, pt);
      setNni(im+ Nx * (j + Ny * km), neighbour[1] * neighbour[5], &nLocal[15], ni, pt);
      setNni(i + Nx * (jp+ Ny * kp), neighbour[2] * neighbour[4], &nLocal[16], ni, pt);
      setNni(i + Nx * (jp+ Ny * km), neighbour[2] * neighbour[5], &nLocal[17], ni, pt);
      setNni(i + Nx * (jm+ Ny * kp), neighbour[3] * neighbour[4], &nLocal[18], ni, pt);
      setNni(i + Nx * (jm+ Ny * km), neighbour[3] * neighbour[5], &nLocal[19], ni, pt);
    }
    
    // Alteração 22: Acesso a neighbourKind aninhado em 'neighbourhood'
    if (params->neighbourhood.neighbourKind == 3) {
      setNni(ip + Nx * (jp + Ny * kp), neighbour[0] * neighbour[2] * neighbour[4], &nLocal[20], ni, pt);
      setNni(ip + Nx * (jp + Ny * km), neighbour[0] * neighbour[2] * neighbour[5], &nLocal[21], ni, pt);
      setNni(ip + Nx * (jm + Ny * kp), neighbour[0] * neighbour[3] * neighbour[4], &nLocal[22], ni, pt);
      setNni(ip + Nx * (jm + Ny * km), neighbour[0] * neighbour[3] * neighbour[5], &nLocal[23], ni, pt);
      setNni(im + Nx * (jp + Ny * kp), neighbour[1] * neighbour[2] * neighbour[4], &nLocal[24], ni, pt);
      setNni(im + Nx * (jp + Ny * km), neighbour[1] * neighbour[2] * neighbour[5], &nLocal[25], ni, pt);
      setNni(im + Nx * (jm + Ny * kp), neighbour[1] * neighbour[3] * neighbour[4], &nLocal[26], ni, pt);
      setNni(im + Nx * (jm + Ny * km), neighbour[1] * neighbour[3] * neighbour[5], &nLocal[27], ni, pt);
    }
    Etot += geometry->latice_Potential(nLocal);
  }
  return Etot / Nt;
}

void EvolveN::Monte_Carlo_Step(float &ang_var, gsl_rng **r) {
  // Alteração 23: Acesso a T aninhado em 'potential'
  float T = params->potential.T;
  
  // Alteração 24: Acesso a dimensões da grade aninhadas em 'lattice'
  static int Nt = params->lattice.Nx * params->lattice.Ny * params->lattice.Nz;
  static int vallid = 0;
  static float *myNi = (float *)calloc(Nt * 3, sizeof(float));
  int acceptance;
  int idx;
  static int time = 0;
  acceptance = 0;
  int iBoxSize = 2;
  int jBoxSize = 2;
  int kBoxSize = 2;
  // Alteração 25, 26, 27: Acesso a dimensões da grade aninhadas em 'lattice'
  int Ni = params->lattice.Nx / 2 + params->lattice.Nx % 2;
  int Nj = params->lattice.Ny / 2 + params->lattice.Ny % 2;
  int Nk = params->lattice.Nz / 2 + params->lattice.Nz % 2;
  int Ntt = Ni * Nj * Nk;
  static bool setUp = true;
  if (setUp) {
    for (int ii = 0; ii < Nt; ii++)
      if (pt[ii])
        vallid++;
    setUp = false;
  }
  for (int i = 0; i < Nt * 3; i++)
    myNi[i] = ni[i];
#pragma omp parallel num_threads(omp_get_max_threads())
  for (int tick = 0; tick < 8; tick++) {
#pragma omp for simd reduction(+ : acceptance) schedule(simd : dynamic)
    for (int nt = 0; nt < Ntt; nt++) {
      float E_new, E_old, rotation_type, va, ranVal;
      int i, j, k;
      float nNew[3];
      int thread = omp_get_thread_num();
      int iBox = nt % Ni;
      int jBox = (nt / Ni) % Nj;
      int kBox = nt / Ni / Nj;
      int di = tick % 2;
      int dj = (tick / 2) % 2;
      int dk = tick / 4;
      
      // Alteração 28, 29: Acesso a neighbourKind aninhado em 'neighbourhood'
      const static int nti = params->neighbourhood.neighbourKind == 2 ? 20 : (params->neighbourhood.neighbourKind == 3 ? +28 : 8);
      nni nLocal[nti];

      // Alteração 30, 31, 32: Acesso a dimensões da grade aninhadas em 'lattice'
      i = di + iBox * iBoxSize;
      if (i >= params->lattice.Nx)
        continue;
      j = dj + jBox * jBoxSize;
      if (j >= params->lattice.Ny)
        continue;
      k = dk + kBox * kBoxSize;
      if (k >= params->lattice.Nz)
        continue;
        
      if (pti(i, j, k) == 0)
        continue;
        
      int im = (i - 1), ip = (i + 1);
      int jm = (j - 1), jp = (j + 1);
      int km = (k - 1), kp = (k + 1);
      
      // Alteração 33-38: Acesso a funções de fronteira e dimensões aninhadas em 'lattice'
      int neighbour[6] = {
          params->lattice.XBound(ip, params->lattice.Nx), params->lattice.XBound(im, params->lattice.Nx),
          params->lattice.YBound(jp, params->lattice.Ny), params->lattice.YBound(jm, params->lattice.Ny),
          params->lattice.ZBound(kp, params->lattice.Nz), params->lattice.ZBound(km, params->lattice.Nz)};
          
      setNni(i + Nx * (j + Ny * k ), 1, &nLocal[0], ni, pt);
      setNni(ip+ Nx * (j + Ny * k ), neighbour[0], &nLocal[1], ni, pt);
      setNni(im+ Nx * (j + Ny * k ), neighbour[1], &nLocal[2], ni, pt);
      setNni(i + Nx * (jp+ Ny * k ), neighbour[2], &nLocal[3], ni, pt);
      setNni(i + Nx * (jm+ Ny * k ), neighbour[3], &nLocal[4], ni, pt);
      setNni(i + Nx * (j + Ny * kp), neighbour[4], &nLocal[5], ni, pt);
      setNni(i + Nx * (j + Ny * km), neighbour[5], &nLocal[6], ni, pt);
      if (nLocal[0].pt > 1) {
        nLocal[7].x = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 0];
        nLocal[7].y = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 1];
        nLocal[7].z = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 2];
        nLocal[7].pt = 1;
        //~ printf("%g %g %g %d\n",nLocal[7].x,nLocal[7].y,nLocal[7].z,nLocal[7].pt );
      }
      
      // Alteração 39: Acesso a neighbourKind aninhado em 'neighbourhood'
      if (params->neighbourhood.neighbourKind > 1) {
        setNni(ip+ Nx * (jp+ Ny * k ), neighbour[0] * neighbour[2], &nLocal[ 8], ni, pt);
        setNni(ip+ Nx * (jm+ Ny * k ), neighbour[0] * neighbour[3], &nLocal[ 9], ni, pt);
        setNni(ip+ Nx * (j + Ny * kp), neighbour[0] * neighbour[4], &nLocal[10], ni, pt);
        setNni(ip+ Nx * (j + Ny * km), neighbour[0] * neighbour[5], &nLocal[11], ni, pt);
        setNni(im+ Nx * (jp+ Ny * k ), neighbour[1] * neighbour[2], &nLocal[12], ni, pt);
        setNni(im+ Nx * (jm+ Ny * k ), neighbour[1] * neighbour[3], &nLocal[13], ni, pt);
        setNni(im+ Nx * (j + Ny * kp), neighbour[1] * neighbour[4], &nLocal[14], ni, pt);
        setNni(im+ Nx * (j + Ny * km), neighbour[1] * neighbour[5], &nLocal[15], ni, pt);
        setNni(i + Nx * (jp+ Ny * kp), neighbour[2] * neighbour[4], &nLocal[16], ni, pt);
        setNni(i + Nx * (jp+ Ny * km), neighbour[2] * neighbour[5], &nLocal[17], ni, pt);
        setNni(i + Nx * (jm+ Ny * kp), neighbour[3] * neighbour[4], &nLocal[18], ni, pt);
        setNni(i + Nx * (jm+ Ny * km), neighbour[3] * neighbour[5], &nLocal[19], ni, pt);
      }
      
      // Alteração 40: Acesso a neighbourKind aninhado em 'neighbourhood'
      if (params->neighbourhood.neighbourKind == 3) {
        setNni(ip + Nx * (jp + Ny * kp), neighbour[0] * neighbour[2] * neighbour[4], &nLocal[20], ni, pt);
        setNni(ip + Nx * (jp + Ny * km), neighbour[0] * neighbour[2] * neighbour[5], &nLocal[21], ni, pt);
        setNni(ip + Nx * (jm + Ny * kp), neighbour[0] * neighbour[3] * neighbour[4], &nLocal[22], ni, pt);
        setNni(ip + Nx * (jm + Ny * km), neighbour[0] * neighbour[3] * neighbour[5], &nLocal[23], ni, pt);
        setNni(im + Nx * (jp + Ny * kp), neighbour[1] * neighbour[2] * neighbour[4], &nLocal[24], ni, pt);
        setNni(im + Nx * (jp + Ny * km), neighbour[1] * neighbour[2] * neighbour[5], &nLocal[25], ni, pt);
        setNni(im + Nx * (jm + Ny * kp), neighbour[1] * neighbour[3] * neighbour[4], &nLocal[26], ni, pt);
        setNni(im + Nx * (jm + Ny * km), neighbour[1] * neighbour[3] * neighbour[5], &nLocal[27], ni, pt);
      }
      //~ va=gsl_rng_uniform(r)*ang_var;
      va = (2 * gsl_rng_uniform(r[thread]) - 1) * ang_var;
      // Create a rotation candidate
      rotation_type = 0.7;//gsl_rng_uniform(r[thread]);
      {
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
      }

      // Old energy calculation
      E_old = geometry->latice_Potential(nLocal);
      // New energy calculation
      nLocal[0].x = nNew[0];
      nLocal[0].y = nNew[1];
      nLocal[0].z = nNew[2];
      E_new = geometry->latice_Potential(nLocal);

      // test new config
      if (gsl_rng_uniform(r[thread]) < exp(-(E_new - E_old) / T)) {
        nix(i, j, k) = nNew[0];
        niy(i, j, k) = nNew[1];
        niz(i, j, k) = nNew[2];
        acceptance++;
      }
    }
  }

      // printf("- %d %g %d %g\n",acceptance, 1.0*acceptance/vallid,vallid,ang_var);
  //~   #pragma omp barrier
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

// Alteração 41: Assinatura da função check_Points (adicionando const&)
void Evolve::check_Points(int *pt, Parameters params) {
  VallidPoints = 0;
  // Alteração 42, 43, 44: Acesso a dimensões da grade aninhadas em 'lattice'
  int Nx = params.lattice.Nx;
  int Ny = params.lattice.Ny;
  int Nz = params.lattice.Nz;
  for (int k = 0; k < Nz; k++) {
    for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
        if (pti(i, j, k))
          VallidPoints++;
      }
    }
  }
}