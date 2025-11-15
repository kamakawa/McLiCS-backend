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

// Função auxiliar inline para configurar vizinhos
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

// Função auxiliar para configurar todos os vizinhos (elimina duplicação)
void setupNeighbors(int i, int j, int k, int Nx, int Ny, int Nz, 
                   Parameters* params, float* ni, int* pt, Geometry* geometry, 
                   nni* nLocal) {
  int im = (i - 1), ip = (i + 1);
  int jm = (j - 1), jp = (j + 1);
  int km = (k - 1), kp = (k + 1);
  
  int neighbour[6] = {
      params->XBound(ip, Nx), params->XBound(im, Nx),
      params->YBound(jp, Ny), params->YBound(jm, Ny),
      params->ZBound(kp, Nz), params->ZBound(km, Nz)};
  
  // Vizinhos primários (6 primeiros)
  setNni(i + Nx * (j + Ny * k ), 1, &nLocal[0], ni, pt);
  setNni(ip + Nx * (j + Ny * k ), neighbour[0], &nLocal[1], ni, pt);
  setNni(im + Nx * (j + Ny * k ), neighbour[1], &nLocal[2], ni, pt);
  setNni(i + Nx * (jp + Ny * k ), neighbour[2], &nLocal[3], ni, pt);
  setNni(i + Nx * (jm + Ny * k ), neighbour[3], &nLocal[4], ni, pt);
  setNni(i + Nx * (j + Ny * kp), neighbour[4], &nLocal[5], ni, pt);
  setNni(i + Nx * (j + Ny * km), neighbour[5], &nLocal[6], ni, pt);
  
  // Superfície (se aplicável)
  if (nLocal[0].pt > 1) {
    nLocal[7].x = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 0];
    nLocal[7].y = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 1];
    nLocal[7].z = geometry->ns[(i + Nx * (j + Ny * k)) * 3 + 2];
    nLocal[7].pt = 1;
  }
  
  // Vizinhos secundários (neighbourKind > 1)
  if (params->neighbourKind > 1) {
    setNni(ip + Nx * (jp + Ny * k ), neighbour[0] * neighbour[2], &nLocal[8], ni, pt);
    setNni(ip + Nx * (jm + Ny * k ), neighbour[0] * neighbour[3], &nLocal[9], ni, pt);
    setNni(ip + Nx * (j + Ny * kp), neighbour[0] * neighbour[4], &nLocal[10], ni, pt);
    setNni(ip + Nx * (j + Ny * km), neighbour[0] * neighbour[5], &nLocal[11], ni, pt);
    setNni(im + Nx * (jp + Ny * k ), neighbour[1] * neighbour[2], &nLocal[12], ni, pt);
    setNni(im + Nx * (jm + Ny * k ), neighbour[1] * neighbour[3], &nLocal[13], ni, pt);
    setNni(im + Nx * (j + Ny * kp), neighbour[1] * neighbour[4], &nLocal[14], ni, pt);
    setNni(im + Nx * (j + Ny * km), neighbour[1] * neighbour[5], &nLocal[15], ni, pt);
    setNni(i + Nx * (jp + Ny * kp), neighbour[2] * neighbour[4], &nLocal[16], ni, pt);
    setNni(i + Nx * (jp + Ny * km), neighbour[2] * neighbour[5], &nLocal[17], ni, pt);
    setNni(i + Nx * (jm + Ny * kp), neighbour[3] * neighbour[4], &nLocal[18], ni, pt);
    setNni(i + Nx * (jm + Ny * km), neighbour[3] * neighbour[5], &nLocal[19], ni, pt);
  }
  
  // Vizinhos terciários (neighbourKind == 3)
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
}

float EvolveN::energy_calculator() {
  static const int Nt = params->Nx * params->Ny * params->Nz;
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
  const int nti = params->neighbourKind == 2 ? 20 : (params->neighbourKind == 3 ? 28 : 8);
  
#pragma omp parallel for reduction(+ : Etot) schedule(dynamic) num_threads(omp_get_max_threads())
  for (int idx = 0; idx < Nt; idx++) {
    nni nLocal[28];

    int i = idx % Nx;
    int j = (idx / Nx) % Ny;
    int k = idx / (Nx * Ny);
    
    // Verificação de limites
    if (i >= Nx || i < 0 || j >= Ny || j < 0 || k >= Nz || k < 0) {
      printf("Tem cachorro nesse mato (%d %d %d)\n", i, j, k);
      fflush(stdout);
      exit(1);
    }

    // Configura todos os vizinhos
    setupNeighbors(i, j, k, Nx, Ny, Nz, params, ni, pt, geometry, nLocal);
    
    // Calcula energia para este ponto
    Etot += geometry->lattice_Potential(nLocal);
  }
  
  return Etot / Nt;
}

void EvolveN::Monte_Carlo_Step(float &ang_var, gsl_rng **r) {
  float T = params->T;
  static int Nt = params->Nx * params->Ny * params->Nz;
  static int valid = 0;
  static bool setUp = true;
  
  if (setUp) {
    for (int ii = 0; ii < Nt; ii++) {
      if (pt[ii])
        valid++;
    }
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
  const int nti = params->neighbourKind == 2 ? 20 : (params->neighbourKind == 3 ? 28 : 8);

#pragma omp parallel num_threads(omp_get_max_threads())
  {
    for (int tick = 0; tick < 8; tick++) {
#pragma omp for reduction(+ : acceptance) schedule(dynamic)
      for (int nt = 0; nt < Ntt; nt++) {
        int thread = omp_get_thread_num();
        
        // Calcula coordenadas
        int iBox = nt % Ni;
        int jBox = (nt / Ni) % Nj;
        int kBox = nt / (Ni * Nj);
        int di = tick % 2;
        int dj = (tick / 2) % 2;
        int dk = tick / 4;
        
        int i = di + iBox * iBoxSize;
        int j = dj + jBox * jBoxSize;
        int k = dk + kBox * kBoxSize;
        
        // Verifica limites e se o ponto é válido
        if (i >= Nx || j >= Ny || k >= Nz || pti(i, j, k) == 0) {
          continue;
        }

        nni nLocal[28];
        setupNeighbors(i, j, k, Nx, Ny, Nz, params, ni, pt, geometry, nLocal);

        // Gera nova orientação
        float va = (2 * gsl_rng_uniform(r[thread]) - 1) * ang_var;
        float rotation_type = 0.7; // CORREÇÃO: Valor hardcoded - manter compatibilidade
        float nNew[3];
        
        // Aplica rotação (mantendo lógica original)
        if (rotation_type < 0.333f) {
          nNew[0] = nix(i, j, k);
          nNew[1] = niy(i, j, k) * cosf(M_PI * va) + niz(i, j, k) * sinf(M_PI * va);
          nNew[2] = niz(i, j, k) * cosf(M_PI * va) - niy(i, j, k) * sinf(M_PI * va);
        } else if (rotation_type < 0.666f) {
          nNew[0] = nix(i, j, k) * cosf(M_PI * va) - niz(i, j, k) * sinf(M_PI * va);
          nNew[1] = niy(i, j, k);
          nNew[2] = niz(i, j, k) * cosf(M_PI * va) + nix(i, j, k) * sinf(M_PI * va);
        } else {
          nNew[0] = nix(i, j, k) * cosf(M_PI * va) + niy(i, j, k) * sinf(M_PI * va);
          nNew[1] = niy(i, j, k) * cosf(M_PI * va) - nix(i, j, k) * sinf(M_PI * va);
          nNew[2] = niz(i, j, k);
        }

        // Calcula energias
        float E_old = geometry->lattice_Potential(nLocal);
        nLocal[0].x = nNew[0];
        nLocal[0].y = nNew[1];
        nLocal[0].z = nNew[2];
        float E_new = geometry->lattice_Potential(nLocal);

        // Teste de Metropolis
        if (gsl_rng_uniform(r[thread]) < expf(-(E_new - E_old) / T)) {
          nix(i, j, k) = nNew[0];
          niy(i, j, k) = nNew[1];
          niz(i, j, k) = nNew[2];
          acceptance++;
        }
      }
    }
  }

  // Ajusta ang_var baseado na taxa de aceitação
  float acceptance_rate = 1.0f * acceptance / valid;
  if (acceptance_rate < 0.5f) {
    ang_var *= 0.99f;
    if (ang_var < 0.01f)
      ang_var = 0.01f;
  } else if (acceptance_rate > 0.5f) {
    ang_var /= 0.99f;
    if (ang_var > 1.0f) {
      ang_var -= 0.5f;
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
        if (pti(i, j, k))
          VallidPoints++;
      }
    }
  }
}