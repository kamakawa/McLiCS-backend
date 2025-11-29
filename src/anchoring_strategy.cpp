#include <math.h>
#include <iostream>
#include <map>

// --- Project Includes ---
#include "../include/anchoring_strategy.h"
#include "../include/parameters.h"

// ========== Implementação do Método Utilitário Compartilhado ==========

void AnchoringStrategy::checkParameter(bool hasValue, const std::string& parameterName, int id, const std::string& strategyName) {
    if (hasValue) {
        std::cout << "Parameter " << parameterName << " is not set. Using standard value.\n";
    } else {
        std::cout << "Parameter " << parameterName << " not defined for the boundary condition "
                  << "#" << id << ".\n";
        std::cout << "The boundary condition " << strategyName << " needs the aforementioned parameter defined.\n";
        std::cout << "Please, set it in your input file, or check if it is mispelled.\n";
        std::cout << "Aborting the program.\n";
        exit(1);
    }
}

// ========== Implementação do Strong Anchoring ==========

void StrongAnchoringStrategy::initialize(Parameters* params, int id) {
    this->params_ = params;
    this->id_ = id;
    
    printf("setting surface %d: %s\n", id, getName().c_str());
    
    try {
        W_ = params->W.at(id);
        std::cout << "W= " << W_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "W", id, getName());
    }

    try {
        phi_s_ = params->phi_s.at(id);
        std::cout << "phi_s= " << phi_s_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "phi_s", id, getName());
    }

    try {
        theta_s_ = params->theta_s.at(id);
        std::cout << "theta_s= " << theta_s_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "theta_s", id, getName());
    }
    printf("\n");
}

float StrongAnchoringStrategy::calculateSurfacePotential(float ni[3], float s[3], Parameters* params) {
    const float toPi = M_PI / 180.0f;
    
    float n_s[3] = {
        cos(toPi * phi_s_) * sin(toPi * theta_s_),
        sin(toPi * phi_s_) * sin(toPi * theta_s_),
        cos(toPi * theta_s_)
    };

    float nij = ni[0] * n_s[0] + ni[1] * n_s[1] + ni[2] * n_s[2];
    
    return -W_ * nij * nij;
}

// ========== Implementação do Strong Anchoring GHRL ==========

void StrongAnchoringGHRLStrategy::initialize(Parameters* params, int id) {
    this->params_ = params;
    this->id_ = id;
    
    printf("setting surface %d: %s\n", id, getName().c_str());
    
    try {
        W_ = params->W.at(id);
        std::cout << "W= " << W_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "W", id, getName());
    }

    try {
        phi_s_ = params->phi_s.at(id);
        std::cout << "phi_s= " << phi_s_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "phi_s", id, getName());
    }

    try {
        theta_s_ = params->theta_s.at(id);
        std::cout << "theta_s= " << theta_s_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "theta_s", id, getName());
    }
    printf("\n");
}

float StrongAnchoringGHRLStrategy::calculateSurfacePotential(float ni[3], float s[3], Parameters* params) {
    const float toPi = M_PI / 180.0f;
    
    float n_s[3] = {
        cos(toPi * phi_s_) * sin(toPi * theta_s_),
        sin(toPi * phi_s_) * sin(toPi * theta_s_),
        cos(toPi * theta_s_)
    };

    float nij = ni[0] * n_s[0] + ni[1] * n_s[1] + ni[2] * n_s[2];
    
    return -W_ * nij * nij;
}

// ========== Implementação do Rapini-Papoular ==========

void RPAnchoringStrategy::initialize(Parameters* params, int id) {
    this->params_ = params;
    this->id_ = id;
    
    printf("setting surface %d: %s\n", id, getName().c_str());
    
    try {
        W_ = params->W.at(id);
        std::cout << "W= " << W_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "W", id, getName());
    }

    try {
        phi_s_ = params->phi_s.at(id);
        std::cout << "phi_s= " << phi_s_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "phi_s", id, getName());
    }

    try {
        theta_s_ = params->theta_s.at(id);
        std::cout << "theta_s= " << theta_s_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "theta_s", id, getName());
    }

    if (params->neighbourKind == 2) {
        W_ *= 4;
        std::cout << "W rescaled by 4 to " << W_ << " for neighbourKind=2\n";
    }
    if (params->neighbourKind == 3) {
        W_ *= 5;
        std::cout << "W rescaled by 5 to " << W_ << " for neighbourKind=3\n";
    }
    printf("\n");
}

float RPAnchoringStrategy::calculateSurfacePotential(float ni[3], float s[3], Parameters* params) {
    const float toPi = M_PI / 180.0f;
    
    float n_s[3] = {
        cos(toPi * phi_s_) * sin(toPi * theta_s_),
        sin(toPi * phi_s_) * sin(toPi * theta_s_),
        cos(toPi * theta_s_)
    };

    float nij = ni[0] * n_s[0] + ni[1] * n_s[1] + ni[2] * n_s[2];
    
    return -W_ * nij * nij;
}

// ========== Implementação do Rapini-Papoular GHRL ==========

void RPAnchoringGHRLStrategy::initialize(Parameters* params, int id) {
    this->params_ = params;
    this->id_ = id;
    
    printf("setting surface %d: %s\n", id, getName().c_str());
    
    try {
        W_ = params->W.at(id);
        std::cout << "W= " << W_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "W", id, getName());
    }

    try {
        phi_s_ = params->phi_s.at(id);
        std::cout << "phi_s= " << phi_s_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "phi_s", id, getName());
    }

    try {
        theta_s_ = params->theta_s.at(id);
        std::cout << "theta_s= " << theta_s_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "theta_s", id, getName());
    }

    if (params->neighbourKind == 2) {
        W_ *= 4;
        std::cout << "W rescaled by 4 to " << W_ << " for neighbourKind=2\n";
    }
    if (params->neighbourKind == 3) {
        W_ *= 5;
        std::cout << "W rescaled by 5 to " << W_ << " for neighbourKind=3\n";
    }
    printf("\n");
}

float RPAnchoringGHRLStrategy::calculateSurfacePotential(float ni[3], float s[3], Parameters* params) {
    const float toPi = M_PI / 180.0f;
    
    float nj[3] = {
        cos(toPi * phi_s_) * sin(toPi * theta_s_),
        sin(toPi * phi_s_) * sin(toPi * theta_s_),
        cos(toPi * theta_s_)
    };

    const float el = params->ghrl_lambda;
    const float em = params->ghrl_mu;
    const float en = params->ghrl_nu;
    const float er = params->ghrl_rho;
    const float es = params->ghrl_sigma;
    
    float v15 = 1.5;
    float v05 = 0.5;

    float ai = ni[0] * s[0] + ni[1] * s[1] + ni[2] * s[2];
    float aj = nj[0] * s[0] + nj[1] * s[1] + nj[2] * s[2];
    float nij = ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2];
    
    float pij = v15 * nij * nij - v05;
    
    float cross = (ni[2] * nj[1] - ni[1] * nj[2]) * s[0] + 
                  (ni[0] * nj[2] - ni[2] * nj[0]) * s[1] + 
                  (ni[1] * nj[0] - ni[0] * nj[1]) * s[2];

    float E1 = ((v15 * ai * ai) + (v15 * aj * aj) - 1);
    
    return W_ * ((E1 * (er * pij + el) + em * (ai * aj * nij) - (1.0 / 9.0)) + 
                en * pij + es * (nij > 0 ? 1 : -1) * cross);
}

// ========== Implementação do Fournier-Galatola ==========

void FGAnchoringStrategy::initialize(Parameters* params, int id) {
    this->params_ = params;
    this->id_ = id;
    
    printf("setting surface %d: %s\n", id, getName().c_str());

    try {
        W_ = params->W.at(id);
        std::cout << "W= " << W_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "W", id, getName());
    }

    if (params->neighbourKind == 2) {
        W_ *= 4;
        std::cout << "W rescaled by 4 to " << W_ << " for neighbourKind=2\n";
    }
    if (params->neighbourKind == 3) {
        W_ *= 5;
        std::cout << "W rescaled by 5 to " << W_ << " for neighbourKind=3\n";
    }
    printf("\n");
}

float FGAnchoringStrategy::calculateSurfacePotential(float ni[3], float s[3], Parameters* params) {
    float nij = ni[0] * s[0] + ni[1] * s[1] + ni[2] * s[2];
    
    return +W_ * nij * nij;
}

// ========== Implementação do Fournier-Galatola GHRL ==========

void FGAnchoringGHRLStrategy::initialize(Parameters* params, int id) {
    this->params_ = params;
    this->id_ = id;
    
    printf("setting surface %d: %s\n", id, getName().c_str());
    
    try {
        W_ = params->W.at(id);
        std::cout << "W= " << W_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "W", id, getName());
    }

    if (params->neighbourKind == 2) {
        W_ *= 4;
        std::cout << "W rescaled by 4 to " << W_ << " for neighbourKind=2\n";
    }
    if (params->neighbourKind == 3) {
        W_ *= 6;
        std::cout << "W rescaled by 6 to " << W_ << " for neighbourKind=3\n";
    }
    printf("\n");
}

float FGAnchoringGHRLStrategy::calculateSurfacePotential(float ni[3], float s[3], Parameters* params) {
    const float el = params->ghrl_lambda;
    const float em = params->ghrl_mu;
    const float en = params->ghrl_nu;
    const float er = params->ghrl_rho;
    const float es = params->ghrl_sigma;
    
    float v15 = 1.5;
    float v05 = 0.5;

    float mod = sqrtf(fabs(ni[0] * s[0] + ni[1] * s[1] + ni[2] * s[2]));
    float nj[3] = {ni[0] - s[0] * mod, ni[1] - s[1] * mod, ni[2] - s[2] * mod};
    
    mod = sqrtf(fabs(nj[0] * nj[0] + nj[1] * nj[1] + nj[2] * nj[2]));
    
    float scale = (mod > 0 ? mod : 1);
    nj[0] /= scale;
    nj[1] /= scale;
    nj[2] /= scale;

    float ai = ni[0] * s[0] + ni[1] * s[1] + ni[2] * s[2];
    float aj = nj[0] * s[0] + nj[1] * s[1] + nj[2] * s[2];
    float nij = ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2];
    
    float pij = v15 * nij * nij - v05;
    
    float cross = (ni[2] * nj[1] - ni[1] * nj[2]) * s[0] + 
                  (ni[0] * nj[2] - ni[2] * nj[0]) * s[1] + 
                  (ni[1] * nj[0] - ni[0] * nj[1]) * s[2];

    float E1 = ((v15 * ai * ai) + (v15 * aj * aj) - 1);
    
    return W_ * ((E1 * (er * pij + el) + em * (ai * aj * nij) - (1.0 / 9.0)) + 
                en * pij + es * (nij > 0 ? 1 : -1) * cross);
}

// ========== Implementação do Homeotropic ==========

void HomeotropicAnchoringStrategy::initialize(Parameters* params, int id) {
    this->params_ = params;
    this->id_ = id;
    
    printf("setting surface %d: %s\n", id, getName().c_str());
    
    try {
        W_ = params->W.at(id);
        std::cout << "W= " << W_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "W", id, getName());
    }

    if (params->neighbourKind == 2) {
        W_ *= 4;
        std::cout << "W rescaled by 4 to " << W_ << " for neighbourKind=2\n";
    }
    if (params->neighbourKind == 3) {
        W_ *= 5;
        std::cout << "W rescaled by 5 to " << W_ << " for neighbourKind=3\n";
    }
    printf("\n");
}

float HomeotropicAnchoringStrategy::calculateSurfacePotential(float ni[3], float s[3], Parameters* params) {
    float nij = ni[0] * s[0] + ni[1] * s[1] + ni[2] * s[2];
    
    return -W_ * nij * nij;
}

// ========== Implementação do Homeotropic GHRL ==========

void HomeotropicAnchoringGHRLStrategy::initialize(Parameters* params, int id) {
    this->params_ = params;
    this->id_ = id;
    
    printf("setting surface %d: %s\n", id, getName().c_str());
    
    try {
        W_ = params->W.at(id);
        std::cout << "W= " << W_ << ".\n";
    } catch (const std::out_of_range &) {
        checkParameter(false, "W", id, getName());
    }

    if (params->neighbourKind == 2) {
        W_ *= 4;
        std::cout << "W rescaled by 4 to " << W_ << " for neighbourKind=2\n";
    }
    printf("\n");
}

float HomeotropicAnchoringGHRLStrategy::calculateSurfacePotential(float ni[3], float s[3], Parameters* params) {
    const float el = params->ghrl_lambda;
    const float em = params->ghrl_mu;
    const float en = params->ghrl_nu;
    const float er = params->ghrl_rho;
    const float es = params->ghrl_sigma;
    
    float v15 = 1.5;
    float v05 = 0.5;

    float ai = ni[0] * s[0] + ni[1] * s[1] + ni[2] * s[2];
    float aj = 1.0;
    float nij = ni[0] * s[0] + ni[1] * s[1] + ni[2] * s[2];
    
    float pij = v15 * nij * nij - v05;
    
    float cross = (ni[2] * s[1] - ni[1] * s[2]) * s[0] + 
                  (ni[0] * s[2] - ni[2] * s[0]) * s[1] + 
                  (ni[1] * s[0] - ni[0] * s[1]) * s[2];

    float E1 = ((v15 * ai * ai) + (v15 * aj * aj) - 1);

    return W_ * ((E1 * (er * pij + el) + em * (ai * aj * nij) - (1.0 / 9.0)) + 
                en * pij + es * (nij > 0 ? 1 : -1) * cross);
}

// ========== Factory Method ==========

AnchoringStrategy* AnchoringStrategyFactory::create(const std::string& anchoringType) {
    if (anchoringType == "strong") {
        return new StrongAnchoringStrategy();
    } else if (anchoringType == "strong_ghrl") {
        return new StrongAnchoringGHRLStrategy();
    } else if (anchoringType == "rp") {
        return new RPAnchoringStrategy();
    } else if (anchoringType == "rp_ghrl") {
        return new RPAnchoringGHRLStrategy();
    } else if (anchoringType == "fg") {
        return new FGAnchoringStrategy();
    } else if (anchoringType == "fg_ghrl") {
        return new FGAnchoringGHRLStrategy();
    } else if (anchoringType == "homeotropic") {
        return new HomeotropicAnchoringStrategy();
    } else if (anchoringType == "homeotropic_ghrl") {
        return new HomeotropicAnchoringGHRLStrategy();
    } else {
        std::cerr << "Unknown anchoring type: " << anchoringType << std::endl;
        return nullptr;
    }
}

void AnchoringStrategyFactory::initializeAnchoringStrategies(Parameters* params, std::vector<AnchoringStrategy*>& strategies) {
    std::string anchoring;
    int numSurfaces = strategies.capacity(); // Assume que o vector tem capacidade pré-alocada
    
    for (int ii = 0; ii < numSurfaces; ii++) {
        try {
            anchoring = params->anchoring_type.at(ii);
        } catch (std::out_of_range &) {
            std::cout << "You must define " << numSurfaces << " boundaries.\n"
                      << "Please review your input file.\nAborting the program.\n\n";
            exit(0);
        }

        AnchoringStrategy* strategy = create(anchoring);
        if (strategy) {
            strategy->initialize(params, ii);
            strategies[ii] = strategy;
        } else {
            printf("%s boundary condition is not defined\n", anchoring.c_str());
            exit(2);
        }
    }
}