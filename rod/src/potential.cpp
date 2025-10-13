#include "../include/potential.h"
#include <math.h>

#include <iostream>

#include "../include/define.h"
#include "../include/parameters.h"

// Abre o bloco de namespace
namespace Potential { 

// Alteração 1: Assinatura com const e função no namespace Potential (REMOVIDO 'Potential::')
float Bulk_Energy_Selinger_Pear(const float ni[3], const float nj[3], const Parameters *params, const float rij[3], int nk) {
  float nij, polar_splay;
  nij = ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2];
  polar_splay = +(1 + 2 * nij + nij * nij) * ((nj[0] - ni[0]) * rij[0] + (nj[1] - ni[1]) * rij[1] + (nj[2] - ni[2]) * rij[2]);
  
  // Alteração 2-4: Acesso a parâmetros de Potencial aninhados em 'potential'
  return (-params->potential.A * nij * nij - params->potential.B1 * nij - (params->potential.C * 0.25) * polar_splay);
}

// Alteração 5: Assinatura com const e função no namespace Potential (REMOVIDO 'Potential::')
float Bulk_Energy_Lebwohl_Lasher(const float ni[3], const float nj[3], const Parameters *params, const float rij[3], int nk) {
  float nij = ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2];
  // Alteração 6: Acesso a parâmetro de Potencial aninhado em 'potential'
  return -params->potential.A * 1.5 * nij * nij;
}

// Alteração 7: Assinatura com const e função no namespace Potential (REMOVIDO 'Potential::')
float Bulk_Energy_GHRL(const float ni[3], const float nj[3], const Parameters *params, const float rij[3], int nk) {
  // Alteração 8-17: Acesso a escalas (neighbourhood) e parâmetros GHRL (potential) aninhados
  const float el = ( nk==2 ? params->neighbourhood.lambdaScale: 1)*params->potential.ghrl_lambda;
  const float em = ( nk==2 ? params->neighbourhood.muScale    : 1)*params->potential.ghrl_mu;
  const float en = ( nk==2 ? params->neighbourhood.nuScale    : 1)*params->potential.ghrl_nu;
  const float er = ( nk==2 ? params->neighbourhood.rhoScale   : 1)*params->potential.ghrl_rho;
  const float es = ( nk==2 ? params->neighbourhood.sigmaScale : 1)*params->potential.ghrl_sigma;//((nk - 1) * params->neighbourScale + 2) / (2 * sqrt(nk)) * params->ghrl_sigma;
  float ai = ni[0] * rij[0] + ni[1] * rij[1] + ni[2] * rij[2];
  float aj = nj[0] * rij[0] + nj[1] * rij[1] + nj[2] * rij[2];
  float nij = ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2];
  float pij = ((float)1.5) * nij * nij - ((float)0.5);
  float cross = (ni[2] * nj[1] - ni[1] * nj[2]) * rij[0] 
              + (ni[0] * nj[2] - ni[2] * nj[0]) * rij[1] 
              + (ni[1] * nj[0] - ni[0] * nj[1]) * rij[2];

  float E1 = ((1.5 * ai * ai) + (1.5 * aj * aj) - 1);
  return ((E1 * (er * pij + el) + em * (ai * aj * nij - 1.0 / 9)) + en * pij + es * (nij > 0 ? 1 : -1) * cross);
}

// Alteração 18: Assinatura com const e função no namespace Potential (REMOVIDO 'Potential::')
float Electric_Potential(const float ni[3], const Parameters *params){
  // Alteração 19-21: Acesso a elecE, elecX, elecY, elecZ aninhados em 'electric'
  float nDotE = params->electric.elecE*(ni[0]*params->electric.elecX+ni[1]*params->electric.elecY+ni[2]*params->electric.elecZ);

  // Alteração 22: Acesso a elecA aninhado em 'electric'
  return -params->electric.elecA*( nDotE*nDotE );
}

// Alteração 23: Função no namespace Potential (REMOVIDO 'Potential::')
int Periodic_Boundary(int &ii, int NN) {
  ii = ((ii + NN) % NN);
  return 1;
}

// Alteração 24: Função no namespace Potential (REMOVIDO 'Potential::')
int Free_Boundary(int &ii, int NN) {
  if (ii == -1 || ii == NN)
    return 0;
  else
    return 1;
}

// Fecha o bloco de namespace
}