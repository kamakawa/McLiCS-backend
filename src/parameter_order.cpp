#include "../include/parameter_order.h"

// --- System Includes ---
#include <gsl/gsl_eigen.h>
#include <math.h>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/parameters.h"

namespace OrderParam {

float Eigen_value_evaluation(float *mat, float *vec) {
  // ============================================================
  // MELHORIA: thread_local para evitar data race caso o usuário
  // paralelize esta rotina no futuro. Mantém performance (cache).
  // ============================================================
  static thread_local gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);
  static thread_local gsl_vector *eval = gsl_vector_alloc(3);
  static thread_local gsl_matrix *evec = gsl_matrix_alloc(3, 3);
  static thread_local gsl_matrix *m = gsl_matrix_alloc(3, 3);

  float temp;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = mat[3 * i + j];
      gsl_matrix_set(m, j, i, temp);
    }
  }

  gsl_eigen_symmv(m, eval, evec, w);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

  for (int i = 0; i < 3; i++) {
    vec[i] = gsl_matrix_get(evec, i, 0);
  }

  return (gsl_vector_get(eval, 0));
}

void Matrice_constructor(float *ni, float *Q, int *pt, Parameters params) {
  // ============================================================
  // BUGFIX CRÍTICO:
  // Nx/Ny/Nz NÃO podem ser static, senão ficam presos ao primeiro
  // params usado. Isso pode gerar resultados errados se params
  // variar entre execuções.
  // ============================================================
  const int Nx = params.Nx;
  const int Ny = params.Ny;
  const int Nz = params.Nz;

  int points = 0;

  for (int i = 0; i < 9; i++) Q[i] = 0.0f;

  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
        if (pti(i, j, k)) {
          Q[0] += 3.0f * (nix(i, j, k) * nix(i, j, k)) - 1.0f;
          Q[1] += 3.0f * (nix(i, j, k) * niy(i, j, k));
          Q[2] += 3.0f * (nix(i, j, k) * niz(i, j, k));
          Q[4] += 3.0f * (niy(i, j, k) * niy(i, j, k)) - 1.0f;
          Q[5] += 3.0f * (niy(i, j, k) * niz(i, j, k));
          points++;
        }
      }
    }
  }

  // ============================================================
  // SAFETY: evita divisão por zero se não houver pontos válidos
  // ============================================================
  if (points == 0) {
    for (int i = 0; i < 9; i++) Q[i] = 0.0f;
    return;
  }

  const float denom = 2.0f * (float)points;

  Q[0] /= denom;
  Q[1] /= denom;
  Q[2] /= denom;
  Q[4] /= denom;
  Q[5] /= denom;

  Q[3] = Q[1];
  Q[6] = Q[2];
  Q[7] = Q[5];
  Q[8] = -Q[0] - Q[4];
}

float lattice_order_parameter(float *ni, int *pt, int i, int j, int k, Parameters params) {
  // ============================================================
  // BUGFIX: Nx/Ny/Nz não podem ser static pelos mesmos motivos.
  // ============================================================
  const int Nx = params.Nx;
  const int Ny = params.Ny;
  const int Nz = params.Nz;

  int points = 0;
  float Q[9];

  for (int ii = 0; ii < 9; ii++) Q[ii] = 0.0f;

  // Loop sobre vizinhos (cubo 3x3x3)
  for (int di = -1; di < 2; di++) {
    for (int dj = -1; dj < 2; dj++) {
      for (int dk = -1; dk < 2; dk++) {

        if (abs(di) + abs(dj) + abs(dk) > 1) continue;

        int ii = i + di;
        int jj = j + dj;
        int kk = k + dk;

        // Verifica condicoes de contorno e validade do ponto
        if ((params.XBound(ii, Nx)) && (params.YBound(jj, Ny)) && (params.ZBound(kk, Nz)) && pti(ii, jj, kk) != 0) {
          int idx_neighbor = (ii + Nx * (jj + Ny * kk)) * 3;

          for (int l = 0; l < 3; l++) {
            for (int m = 0; m <= l; m++) {
              Q[3 * m + l] += 3.0f * (ni[idx_neighbor + l] * ni[idx_neighbor + m]);
            }
            Q[4 * l] -= 1.0f;
          }
          points++;
        }
      }
    }
  }

  // ============================================================
  // SAFETY: se não houver vizinhos válidos, evita divisão por zero
  // ============================================================
  if (points == 0) return 0.0f;

  const float denom = 2.0f * (float)points;

  // Normalizacao
  Q[0] /= denom;
  Q[1] /= denom;
  Q[2] /= denom;
  Q[4] /= denom;
  Q[5] /= denom;
  Q[8] /= denom;

  // Simetrizacao
  Q[3] = Q[1];
  Q[6] = Q[2];
  Q[7] = Q[5];

  // Diagonalizacao
  static thread_local gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);
  static thread_local gsl_vector *eval = gsl_vector_alloc(3);
  static thread_local gsl_matrix *evec = gsl_matrix_alloc(3, 3);
  static thread_local gsl_matrix *m = gsl_matrix_alloc(3, 3);

  for (int l = 0; l < 3; l++) {
    for (int c = 0; c < 3; c++) {
      gsl_matrix_set(m, c, l, Q[3 * l + c]);
    }
  }

  gsl_eigen_symmv(m, eval, evec, w);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

  return (gsl_vector_get(eval, 0));
}

float V_Order_parameter_evaluation(float *mat_b, float *mat_c, float *vec_b, float *vec_c, float PoB) {
  return ((-VMV(mat_b, vec_c) - VMV(mat_c, vec_b) + VMV(mat_c, vec_c) + PoB) / 3.0f);
}

void C_Vector_evaluation(float *vec_n, float *vec_b, float *vec_c) {
  vec_c[0] = vec_n[1] * vec_b[2] - vec_b[1] * vec_n[2];
  vec_c[1] = vec_n[2] * vec_b[0] - vec_b[2] * vec_n[0];
  vec_c[2] = vec_n[0] * vec_b[1] - vec_b[0] * vec_n[1];

  float unit = sqrtf(vec_c[0] * vec_c[0] + vec_c[1] * vec_c[1] + vec_c[2] * vec_c[2]);
  // SAFETY: evita NaN se unit == 0 (caso raro, mas possível)
  if (unit == 0.0f) {
    vec_c[0] = 0.0f; vec_c[1] = 0.0f; vec_c[2] = 0.0f;
    return;
  }
  vec_c[0] /= unit;
  vec_c[1] /= unit;
  vec_c[2] /= unit;
}

float VMV(float *M, float *V) {
  int a, b;
  float prod = 0.0f;

  for (a = 0; a < 3; a++) {
    prod += V[a] * V[a] * M[4 * a];
    for (b = 0; b < a; b++) {
      prod += 2.0f * V[a] * V[b] * M[a + 3 * b];
    }
  }
  return (prod);
}

float Polarization(float *bi, Parameters params) {
  float P_temp[3], modulo;
  P_temp[0] = 0.0f;
  P_temp[1] = 0.0f;
  P_temp[2] = 0.0f;

  const int Nt = params.Nx * params.Ny * params.Nz;

  // SAFETY: evita divisão por zero em casos degenerados
  if (Nt <= 0) return 0.0f;

  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < Nt; i++) {
      P_temp[j] += bi[3 * i + j];
    }
    P_temp[j] /= (float)Nt;
  }

  modulo = sqrtf(P_temp[0] * P_temp[0] + P_temp[1] * P_temp[1] + P_temp[2] * P_temp[2]);
  return (modulo);
}

} // namespace OrderParam