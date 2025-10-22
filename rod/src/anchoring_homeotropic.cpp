//alterei
//Praticamente todo o código foi reescrito para evitar repetição de código, criando funções auxiliares encapsuladas em um namespace anônimo.
//mesmas alterações foram feitas em anchoring_fg.cpp 
//alterações na estrutura, mas a lógica original foi mantida, ou seja, o resultado final não deve ser alterado.
#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <map>
#include <vector>

#include "../include/anchoring.h"
#include "../include/parameters.h"
#include "../include/potential.h"

namespace {

    float getAnchoringW(Parameters* params, int id) {
      try {
          return params->surface.W.at(id);
      } catch (std::out_of_range dummy_var) {
          std::cerr << "Erro: Parâmetro W não encontrado para surface id " << id << std::endl;
          throw std::runtime_error("Parâmetro W de anchoring não configurado");
          // A chamada a check_parameter foi ELIMINADA.
      }
  }
    
    void applyNeighbourScalingHomeotropic(Parameters* params, float& W) {
        if (params->neighbourhood.neighbourKind == 2) {
            W *= 4;
            std::cout << "W reescaled by 4 to " << W << " to acomodate the extra neighbours.\n";
        }
        if (params->neighbourhood.neighbourKind == 3) {
            W *= 5;
            std::cout << "W reescaled by 5 to " << W << " to acomodate the extra neighbours.\n";
        }
    }
    
    void applyNeighbourScalingHomeotropicGHRL(Parameters* params, float& W) {
        if (params->neighbourhood.neighbourKind == 2) {
            W *= 4;
            std::cout << "W reescaled by 4 to " << W << " to acomodate the extra neighbours.\n";
        }
    }

    struct GHRL_Constants {
        float el, em, en, er, es;
    };

    GHRL_Constants getGHRLParams(Parameters* params) {
        return {
            params->potential.ghrl_lambda,
            params->potential.ghrl_mu,
            params->potential.ghrl_nu,
            params->potential.ghrl_rho,
            params->potential.ghrl_sigma
        };
    }
}

Homeotropic_Anchoring::Homeotropic_Anchoring(Parameters *params, int id) {
    this->id = id;
    this->params = params;
    
    printf("seting surface %d: %s\n", id, name);
    this->W = getAnchoringW(params, id);
    std::cout << "W= " << W << ".\n";
    
    applyNeighbourScalingHomeotropic(params, this->W);
    printf("\n");
}

float Homeotropic_Anchoring::surface_potential(float ni[3], float s[3]) {
    float nij = ni[0] * s[0] + ni[1] * s[1] + ni[2] * s[2];
    return -W * nij * nij;
}

Homeotropic_Anchoring_GHRL::Homeotropic_Anchoring_GHRL(Parameters *params, int id) {
    this->id = id;
    this->params = params;
    
    printf("seting surface %d: %s\n", id, name);
    this->W = getAnchoringW(params, id);
    std::cout << "W= " << W << ".\n";
    
    applyNeighbourScalingHomeotropicGHRL(params, this->W);
    printf("\n");
}

float Homeotropic_Anchoring_GHRL::surface_potential(float ni[3], float s[3]) {
    const GHRL_Constants p = getGHRLParams(params);
    
    float v15 = 1.5;
    float v05 = 0.5;

    float ai = ni[0] * s[0] + ni[1] * s[1] + ni[2] * s[2];
    float aj = 1;
    float nij = ni[0] * s[0] + ni[1] * s[1] + ni[2] * s[2];
    float pij = v15 * nij * nij - v05;
    float cross = (ni[2] * s[1] - ni[1] * s[2]) * s[0] + (ni[0] * s[2] - ni[2] * s[0]) * s[1] + (ni[1] * s[0] - ni[0] * s[1]) * s[2];

    float E1 = ((v15 * ai * ai) + (v15 * aj * aj) - 1);
    return W * ((E1 * (p.er * pij + p.el) + p.em * (ai * aj * nij) - (1 / 9)) + p.en * pij + p.es * (nij > 0 ? 1 : -1) * cross);
}
