#include "../include/parameter_order.h"

#include <gsl/gsl_eigen.h>
#include <math.h>
#include <iostream>

#include "../include/define.h"
#include "../include/parameters.h"

// Allocate a fresh workspace per call — avoids sharing state across threads.
// The workspace is tiny (3x3) so allocation cost is negligible.
float Eigen_value_evaluation(float* mat, float* vec) {
  gsl_eigen_symmv_workspace* w    = gsl_eigen_symmv_alloc(3);
  gsl_vector*                eval = gsl_vector_alloc(3);
  gsl_matrix*                evec = gsl_matrix_alloc(3, 3);
  gsl_matrix*                m    = gsl_matrix_alloc(3, 3);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      gsl_matrix_set(m, j, i, mat[3*i + j]);

  gsl_eigen_symmv(m, eval, evec, w);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);
  for (int i = 0; i < 3; i++)
    vec[i] = static_cast<float>(gsl_matrix_get(evec, i, 0));
  const float result = static_cast<float>(gsl_vector_get(eval, 0));

  gsl_matrix_free(m);
  gsl_matrix_free(evec);
  gsl_vector_free(eval);
  gsl_eigen_symmv_free(w);
  return result;
}

void Matrice_constructor(float* ni, float* Q, int* pt, Parameters params) {
  // Use local variables — never static — so a restart with different
  // dimensions picks up the correct values.
  const int Nx = params.Nx;
  const int Ny = params.Ny;
  const int Nz = params.Nz;
  int points = 0;

  for (int i = 0; i < 9; i++) Q[i] = 0;

  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
        if (pti(i, j, k)) {
          Q[0] += 3.f * (nix(i,j,k) * nix(i,j,k)) - 1;
          Q[1] += 3.f * (nix(i,j,k) * niy(i,j,k));
          Q[2] += 3.f * (nix(i,j,k) * niz(i,j,k));
          Q[4] += 3.f * (niy(i,j,k) * niy(i,j,k)) - 1;
          Q[5] += 3.f * (niy(i,j,k) * niz(i,j,k));
          points++;
        }
      }
    }
  }

  const float inv = 1.f / (2.f * points);
  Q[0] *= inv;  Q[1] *= inv;  Q[2] *= inv;
  Q[4] *= inv;  Q[5] *= inv;
  Q[3] = Q[1];  Q[6] = Q[2];  Q[7] = Q[5];
  Q[8] = -Q[0] - Q[4];
}

float lattice_order_parameter(float* ni, int* pt, int i, int j, int k,
                               Parameters params) {
  // Local variables — not static.
  const int Nx = params.Nx;
  const int Ny = params.Ny;
  const int Nz = params.Nz;
  int   points = 0;
  float Q[9]   = {};

  for (int di = -1; di < 2; di++) {
    for (int dj = -1; dj < 2; dj++) {
      for (int dk = -1; dk < 2; dk++) {
        if (abs(di) + abs(dj) + abs(dk) > 1) continue;
        int ii = i + di, jj = j + dj, kk = k + dk;
        if (params.XBound(ii, Nx) && params.YBound(jj, Ny) &&
            params.ZBound(kk, Nz) && pti(ii, jj, kk)) {
          for (int l = 0; l < 3; l++) {
            for (int m2 = 0; m2 <= l; m2++)
              Q[3*m2 + l] += 3.f * (ni[(ii + Nx*(jj + Ny*kk))*3 + l] *
                                     ni[(ii + Nx*(jj + Ny*kk))*3 + m2]);
            Q[4*l] -= 1;
          }
          points++;
        }
      }
    }
  }

  const float inv = 1.f / (2.f * points);
  Q[0] *= inv;  Q[1] *= inv;  Q[2] *= inv;
  Q[4] *= inv;  Q[5] *= inv;  Q[8] *= inv;
  Q[3] = Q[1];  Q[6] = Q[2];  Q[7] = Q[5];

  // Fresh workspace per call — safe for future parallel use.
  gsl_eigen_symmv_workspace* w    = gsl_eigen_symmv_alloc(3);
  gsl_vector*                eval = gsl_vector_alloc(3);
  gsl_matrix*                evec = gsl_matrix_alloc(3, 3);
  gsl_matrix*                m    = gsl_matrix_alloc(3, 3);

  for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
      gsl_matrix_set(m, jj, ii, Q[3*ii + jj]);

  gsl_eigen_symmv(m, eval, evec, w);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);
  const float result = static_cast<float>(gsl_vector_get(eval, 0));

  gsl_matrix_free(m);
  gsl_matrix_free(evec);
  gsl_vector_free(eval);
  gsl_eigen_symmv_free(w);
  return result;
}

float V_Order_parameter_evaluation(float* mat_b, float* mat_c,
                                    float* vec_b, float* vec_c, float PoB) {
  return (-VMV(mat_b, vec_c) - VMV(mat_c, vec_b) +
           VMV(mat_c, vec_c) + PoB) / 3.f;
}

void C_Vector_evaluation(float* vec_n, float* vec_b, float* vec_c) {
  vec_c[0] = vec_n[1]*vec_b[2] - vec_b[1]*vec_n[2];
  vec_c[1] = vec_n[2]*vec_b[0] - vec_b[2]*vec_n[0];
  vec_c[2] = vec_n[0]*vec_b[1] - vec_b[0]*vec_n[1];
  const float unit = sqrtf(vec_c[0]*vec_c[0] + vec_c[1]*vec_c[1] + vec_c[2]*vec_c[2]);
  vec_c[0] /= unit;  vec_c[1] /= unit;  vec_c[2] /= unit;
}

float VMV(float* M, float* V) {
  float prod = 0;
  for (int a = 0; a < 3; a++) {
    prod += V[a] * V[a] * M[4*a];
    for (int b = 0; b < a; b++)
      prod += 2 * V[a] * V[b] * M[a + 3*b];
  }
  return prod;
}

float Polarization(float* bi, Parameters params) {
  float P[3] = {};
  const int Nt = params.Nx * params.Ny * params.Nz;
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < Nt; i++) P[j] += bi[3*i + j];
    P[j] /= Nt;
  }
  return sqrtf(P[0]*P[0] + P[1]*P[1] + P[2]*P[2]);
}
