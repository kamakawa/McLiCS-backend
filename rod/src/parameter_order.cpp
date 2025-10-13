#include "../include/parameter_order.h"

#include <gsl/gsl_eigen.h>
#include <math.h>

#include <iostream>

#include "../include/define.h"
#include "../include/parameters.h"

// Alteração 1: Eigen_value_evaluation movida para o namespace OrderParameters
float OrderParameters::Eigen_value_evaluation(float *mat, float *vec) {
  static gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);
  static gsl_vector *eval = gsl_vector_alloc(3);
  static gsl_matrix *evec = gsl_matrix_alloc(3, 3);
  static gsl_matrix *m = gsl_matrix_alloc(3, 3);
  float temp;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = mat[3 * i + j];
      gsl_matrix_set(m, j, i, temp);
    }
  }

  gsl_eigen_symmv(m, eval, evec, w);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);
  for (int i = 0; i < 3; i++) vec[i] = gsl_matrix_get(evec, i, 0);
  return (gsl_vector_get(eval, 0));
}

// Alteração 2: Matrice_constructor movida para o namespace OrderParameters
void OrderParameters::Matrice_constructor(float *ni, float *Q, int *pt, Parameters params) {
  // Alteração 3, 4, 5: Acesso a dimensões da grade aninhadas em 'lattice'
  static const int Nx = params.lattice.Nx;
  static const int Ny = params.lattice.Ny;
  static const int Nz = params.lattice.Nz;
  int points = 0;
  for (int i = 0; i < 9; i++) Q[i] = 0;
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
        if (pti(i, j, k)) {
          Q[0] += 3.0 * (nix(i, j, k) * nix(i, j, k)) - 1;
          Q[1] += 3.0 * (nix(i, j, k) * niy(i, j, k));
          Q[2] += 3.0 * (nix(i, j, k) * niz(i, j, k));
          Q[4] += 3.0 * (niy(i, j, k) * niy(i, j, k)) - 1;
          Q[5] += 3.0 * (niy(i, j, k) * niz(i, j, k));
          points++;
        }
      }
    }
  }

  Q[0] /= 2.0 * points;
  Q[1] /= 2.0 * points;
  Q[2] /= 2.0 * points;
  Q[4] /= 2.0 * points;
  Q[5] /= 2.0 * points;
  Q[3] = Q[1];
  Q[6] = Q[2];
  Q[7] = Q[5];
  Q[8] = -Q[0] - Q[4];
}

// Alteração 6: lattice_order_parameter movida para o namespace OrderParameters
float OrderParameters::lattice_order_parameter(float *ni, int *pt, int i, int j, int k, Parameters params) {
  // Alteração 7, 8, 9: Acesso a dimensões da grade aninhadas em 'lattice'
  static int Nx = params.lattice.Nx;
  static int Ny = params.lattice.Ny;
  static int Nz = params.lattice.Nz;
  int points = 0;
  float Q[9];
  for (int i = 0; i < 5; i++) Q[i] = 0;

  for (int di = -1; di < 2; di++) {
    for (int dj = -1; dj < 2; dj++) {
      for (int dk = -1; dk < 2; dk++) {
        if (abs(di) + abs(dj) + abs(dk) > 1) continue;
        int ii = i + di;
        int jj = j + dj;
        int kk = k + dk;
        // Alteração 10, 11, 12: Acesso a funções de fronteira e dimensões aninhadas em 'lattice'
        if ((params.lattice.XBound(ii, Nx)) && (params.lattice.YBound(jj, Ny)) && (params.lattice.ZBound(kk, Nz)) && pti(ii, jj, kk) != 0) {
          for (int l = 0; l < 3; l++) {
            for (int m = 0; m <= l; m++) {
              Q[3 * m + l] += 3.0 * (ni[(ii + Nx * (jj + Ny * kk)) * 3 + l] * ni[(ii + Nx * (jj + Ny * kk)) * 3 + m]);
            }
            Q[4 * l] -= 1;
          }
          points++;
        }
      }
    }
  }

  Q[0] /= 2.0 * points;
  Q[1] /= 2.0 * points;
  Q[2] /= 2.0 * points;
  Q[4] /= 2.0 * points;
  Q[5] /= 2.0 * points;
  Q[8] /= 2.0 * points;
  Q[3] = Q[1];
  Q[6] = Q[2];
  Q[7] = Q[5];

  static gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);
  static gsl_vector *eval = gsl_vector_alloc(3);
  static gsl_matrix *evec = gsl_matrix_alloc(3, 3);
  static gsl_matrix *m = gsl_matrix_alloc(3, 3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gsl_matrix_set(m, j, i, Q[3 * i + j]);
    }
  }

  gsl_eigen_symmv(m, eval, evec, w);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);
  return (gsl_vector_get(eval, 0));
}

// Alteração 13: V_Order_parameter_evaluation movida para o namespace OrderParameters
float OrderParameters::V_Order_parameter_evaluation(float *mat_b, float *mat_c, float *vec_b, float *vec_c, float PoB) {
  return ((-VMV(mat_b, vec_c) - VMV(mat_c, vec_b) + VMV(mat_c, vec_c) + PoB) / 3.0);
}

// Alteração 14: C_Vector_evaluation movida para o namespace OrderParameters
void OrderParameters::C_Vector_evaluation(float *vec_n, float *vec_b, float *vec_c) {
  vec_c[0] = vec_n[1] * vec_b[2] - vec_b[1] * vec_n[2];
  vec_c[1] = vec_n[2] * vec_b[0] - vec_b[2] * vec_n[0];
  vec_c[2] = vec_n[0] * vec_b[1] - vec_b[0] * vec_n[1];

  float unit = sqrt(vec_c[0] * vec_c[0] + vec_c[1] * vec_c[1] + vec_c[2] * vec_c[2]);

  vec_c[0] /= unit;
  vec_c[1] /= unit;
  vec_c[2] /= unit;
}

// Alteração 15: VMV movida para o namespace OrderParameters
float OrderParameters::VMV(float *M, float *V) {
  int a, b;
  float prod = 0;

  for (a = 0; a < 3; a++) {
    prod += V[a] * V[a] * M[4 * a];
    for (b = 0; b < a; b++) prod += 2 * V[a] * V[b] * M[a + 3 * b];
  }
  return (prod);
}

// Alteração 16: Polarization movida para o namespace OrderParameters
float OrderParameters::Polarization(float *bi, Parameters params) {
  float P_temp[3], modulo;
  P_temp[0] = 0;
  P_temp[1] = 0;
  P_temp[2] = 0;
  // Alteração 17, 18, 19: Acesso a dimensões da grade aninhadas em 'lattice'
  float Nt = params.lattice.Nx * params.lattice.Ny * params.lattice.Nz;
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < Nt; i++) P_temp[j] += bi[3 * i + j];
    P_temp[j] /= Nt;
  }

  modulo = sqrt(P_temp[0] * P_temp[0] + P_temp[1] * P_temp[1] + P_temp[2] * P_temp[2]);
  return (modulo);
}