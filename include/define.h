#ifndef DEFINE_H_
#define DEFINE_H_

#include "../include/parameters.h"

// ============================================================
// MELHORIA: Centralizar o cálculo do índice e proteger os macros
// com parênteses nos argumentos. Não muda resultado, só torna
// o código mais robusto e reduz repetição.
// ============================================================
#define IDX(i, j, k) ((i) + (Nx) * ((j) + (Ny) * (k)))

#define njx(i, j, k) nj[(IDX((i), (j), (k))) * 3 + 0]
#define njy(i, j, k) nj[(IDX((i), (j), (k))) * 3 + 1]
#define njz(i, j, k) nj[(IDX((i), (j), (k))) * 3 + 2]

#define pti(i, j, k) pt[(IDX((i), (j), (k)))]

#define nix(i, j, k) ni[(IDX((i), (j), (k))) * 3 + 0]
#define niy(i, j, k) ni[(IDX((i), (j), (k))) * 3 + 1]
#define niz(i, j, k) ni[(IDX((i), (j), (k))) * 3 + 2]

#define bix(i, j, k) bi[(IDX((i), (j), (k))) * 3 + 0]
#define biy(i, j, k) bi[(IDX((i), (j), (k))) * 3 + 1]
#define biz(i, j, k) bi[(IDX((i), (j), (k))) * 3 + 2]

#define cix(i, j, k) ci[(IDX((i), (j), (k))) * 3 + 0]
#define ciy(i, j, k) ci[(IDX((i), (j), (k))) * 3 + 1]
#define ciz(i, j, k) ci[(IDX((i), (j), (k))) * 3 + 2]

typedef float (*p_Surface_Potential)(float ni[3], float s[3], Parameters *params, float rij[3]);

struct nni {
  float x;
  float y;
  float z;
  int pt;
};

#endif