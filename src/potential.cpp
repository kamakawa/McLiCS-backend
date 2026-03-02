#include "../include/potential.h"

// --- System Includes ---
#include <math.h>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/parameters.h"

namespace Potential {

// Potencial Selinger-Pear
HD float Bulk_Energy_Selinger_Pear(float ni[3], float nj[3], Parameters *params, float rij[3], int nk) {
  (void)nk;

  float nij = ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2];

  // Mantém exatamente a mesma expressão
  float polar_splay =
      +(1.0f + 2.0f * nij + nij * nij) *
      ((nj[0] - ni[0]) * rij[0] +
       (nj[1] - ni[1]) * rij[1] +
       (nj[2] - ni[2]) * rij[2]);

  return (-params->A * nij * nij - params->B1 * nij - (params->C * 0.25f) * polar_splay);
}

// Potencial Lebwohl-Lasher (Nematico Classico)
HD float Bulk_Energy_Lebwohl_Lasher(float ni[3], float nj[3], Parameters *params, float rij[3], int nk) {
  (void)rij;
  (void)nk;

  float nij = ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2];

  // -A * 1.5 * (n·n)^2
  return -params->A * 1.5f * nij * nij;
}

// Potencial GHRL (Generalizado com escalas para vizinhos)
HD float Bulk_Energy_GHRL(float ni[3], float nj[3], Parameters *params, float rij[3], int nk) {
  // Selecao de escalas baseada na camada de vizinhos (nk)
  const float scale = (nk == 2 ? 1.0f : 1.0f); // mantém estrutura sem mudar nada (nk==2 tem escalas abaixo)

  const float el = (nk == 2 ? params->lambdaScale : 1.0f) * params->ghrl_lambda;
  const float em = (nk == 2 ? params->muScale     : 1.0f) * params->ghrl_mu;
  const float en = (nk == 2 ? params->nuScale     : 1.0f) * params->ghrl_nu;
  const float er = (nk == 2 ? params->rhoScale    : 1.0f) * params->ghrl_rho;
  const float es = (nk == 2 ? params->sigmaScale  : 1.0f) * params->ghrl_sigma;
  (void)scale;

  const float ai  = ni[0] * rij[0] + ni[1] * rij[1] + ni[2] * rij[2];
  const float aj  = nj[0] * rij[0] + nj[1] * rij[1] + nj[2] * rij[2];
  const float nij = ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2];

  constexpr float v15 = 1.5f;
  constexpr float v05 = 0.5f;

  const float pij = v15 * nij * nij - v05;

  const float cross =
      (ni[2] * nj[1] - ni[1] * nj[2]) * rij[0] +
      (ni[0] * nj[2] - ni[2] * nj[0]) * rij[1] +
      (ni[1] * nj[0] - ni[0] * nj[1]) * rij[2];

  const float E1 = (v15 * ai * ai) + (v15 * aj * aj) - 1.0f;

  // Mantém exatamente a mesma expressão do seu código
  return ((E1 * (er * pij + el) + em * (ai * aj * nij - 1.0f / 9.0f)) +
          en * pij + es * (nij > 0.0f ? 1.0f : -1.0f) * cross);
}

// Potencial Eletrico (Interacao dipolo induzido - campo)
HD float Electric_Potential(float ni[3], Parameters *params) {
  const float nDotE =
      params->elecE *
      (ni[0] * params->elecX + ni[1] * params->elecY + ni[2] * params->elecZ);

  return -params->elecA * (nDotE * nDotE);
}

HD int Periodic_Boundary(int &ii, int NN) {
  ii = ((ii + NN) % NN);
  return 1;
}

HD int Free_Boundary(int &ii, int NN) {
  // Retorna 0 se o indice sair da caixa (vazio/vacuo)
  if (ii == -1 || ii == NN)
    return 0;
  else
    return 1;
}

} // namespace Potential