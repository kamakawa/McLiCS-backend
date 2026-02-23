#include "../include/geometry_strategy.h"

// --- System Includes ---
#include <math.h>
#include <string.h>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>   // ← MELHORIA: calloc/free
#include <cstdio>    // ← MELHORIA: fprintf, stderr

// --- Project Includes ---
#include "../include/anchoring_strategy.h"
#include "../include/parameters.h"
#include "../include/potential.h"

// ========== Implementação dos Métodos Compartilhados ==========

// ============================================================
// CORREÇÃO: Assinatura agora usa ponteiro (fullni não tem tamanho fixo).
// Isso é consistente com o uso real (nLocal[8/20/28] no MonteCarlo).
// ============================================================
float GeometryStrategy::calculateNewmanNeighbours(const nni* fullni, Parameters* params) {
    float rij[3];
    double E = 0;
    float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};

    if (!bulk_potential) {
        std::cerr << "Error: bulk_potential not set in GeometryStrategy!" << std::endl;
        return 0.0f;
    }

    // --- Eixo X ---
    if (fullni[1].pt) {
        float nj[3] = {fullni[1].x, fullni[1].y, fullni[1].z};
        rij[0] = 1; rij[1] = 0; rij[2] = 0;
        E += bulk_potential(ni, nj, params, rij, 1);
    }
    if (fullni[2].pt) {
        float nj[3] = {fullni[2].x, fullni[2].y, fullni[2].z};
        rij[0] = -1; rij[1] = 0; rij[2] = 0;
        E += bulk_potential(ni, nj, params, rij, 1);
    }

    // --- Eixo Y ---
    if (fullni[3].pt) {
        float nj[3] = {fullni[3].x, fullni[3].y, fullni[3].z};
        rij[0] = 0; rij[1] = 1; rij[2] = 0;
        E += bulk_potential(ni, nj, params, rij, 1);
    }
    if (fullni[4].pt) {
        float nj[3] = {fullni[4].x, fullni[4].y, fullni[4].z};
        rij[0] = 0; rij[1] = -1; rij[2] = 0;
        E += bulk_potential(ni, nj, params, rij, 1);
    }

    // --- Eixo Z ---
    if (fullni[5].pt) {
        float nj[3] = {fullni[5].x, fullni[5].y, fullni[5].z};
        rij[0] = 0; rij[1] = 0; rij[2] = 1;
        E += bulk_potential(ni, nj, params, rij, 1);
    }
    if (fullni[6].pt) {
        float nj[3] = {fullni[6].x, fullni[6].y, fullni[6].z};
        rij[0] = 0; rij[1] = 0; rij[2] = -1;
        E += bulk_potential(ni, nj, params, rij, 1);
    }
    return (float)E;
}

float GeometryStrategy::calculateSecondNeighbours(const nni* fullni, Parameters* params) {
    float rij[3];
    double E = 0;
    float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
    const float isqrt2 = 0.707106781f;

    if (!bulk_potential) return 0.0f;

    // --- Plano XY ---
    if (fullni[8].pt) {
        float nj[3] = {fullni[8].x, fullni[8].y, fullni[8].z};
        rij[0] = isqrt2; rij[1] = isqrt2; rij[2] = 0;
        E += bulk_potential(ni, nj, params, rij, 2);
    }
    if (fullni[9].pt) {
        float nj[3] = {fullni[9].x, fullni[9].y, fullni[9].z};
        rij[0] = isqrt2; rij[1] = -isqrt2; rij[2] = 0;
        E += bulk_potential(ni, nj, params, rij, 2);
    }
    if (fullni[12].pt) {
        float nj[3] = {fullni[12].x, fullni[12].y, fullni[12].z};
        rij[0] = -isqrt2; rij[1] = isqrt2; rij[2] = 0;
        E += bulk_potential(ni, nj, params, rij, 2);
    }
    if (fullni[13].pt) {
        float nj[3] = {fullni[13].x, fullni[13].y, fullni[13].z};
        rij[0] = -isqrt2; rij[1] = -isqrt2; rij[2] = 0;
        E += bulk_potential(ni, nj, params, rij, 2);
    }

    // --- Plano XZ ---
    if (fullni[10].pt) {
        float nj[3] = {fullni[10].x, fullni[10].y, fullni[10].z};
        rij[0] = isqrt2; rij[1] = 0; rij[2] = isqrt2;
        E += bulk_potential(ni, nj, params, rij, 2);
    }
    if (fullni[11].pt) {
        float nj[3] = {fullni[11].x, fullni[11].y, fullni[11].z};
        rij[0] = isqrt2; rij[1] = 0; rij[2] = -isqrt2;
        E += bulk_potential(ni, nj, params, rij, 2);
    }
    if (fullni[14].pt) {
        float nj[3] = {fullni[14].x, fullni[14].y, fullni[14].z};
        rij[0] = -isqrt2; rij[1] = 0; rij[2] = isqrt2;
        E += bulk_potential(ni, nj, params, rij, 2);
    }
    if (fullni[15].pt) {
        float nj[3] = {fullni[15].x, fullni[15].y, fullni[15].z};
        rij[0] = -isqrt2; rij[1] = 0; rij[2] = -isqrt2;
        E += bulk_potential(ni, nj, params, rij, 2);
    }

    // --- Plano YZ ---
    if (fullni[16].pt) {
        float nj[3] = {fullni[16].x, fullni[16].y, fullni[16].z};
        rij[0] = 0; rij[1] = isqrt2; rij[2] = isqrt2;
        E += bulk_potential(ni, nj, params, rij, 2);
    }
    if (fullni[17].pt) {
        float nj[3] = {fullni[17].x, fullni[17].y, fullni[17].z};
        rij[0] = 0; rij[1] = isqrt2; rij[2] = -isqrt2;
        E += bulk_potential(ni, nj, params, rij, 2);
    }
    if (fullni[18].pt) {
        float nj[3] = {fullni[18].x, fullni[18].y, fullni[18].z};
        rij[0] = 0; rij[1] = -isqrt2; rij[2] = isqrt2;
        E += bulk_potential(ni, nj, params, rij, 2);
    }
    if (fullni[19].pt) {
        float nj[3] = {fullni[19].x, fullni[19].y, fullni[19].z};
        rij[0] = 0; rij[1] = -isqrt2; rij[2] = -isqrt2;
        E += bulk_potential(ni, nj, params, rij, 2);
    }
    return (float)E;
}

float GeometryStrategy::calculateThirdNeighbours(const nni* fullni, Parameters* params) {
    float rij[3];
    double E = 0;
    float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
    const float isqrt3 = 0.577350269f;

    if (!bulk_potential) return 0.0f;

    // --- Vertices Superiores (Z positivo) ---
    if (fullni[20].pt) {
        float nj[3] = {fullni[20].x, fullni[20].y, fullni[20].z};
        rij[0] = isqrt3; rij[1] = isqrt3; rij[2] = isqrt3;
        E += bulk_potential(ni, nj, params, rij, 3);
    }
    if (fullni[22].pt) {
        float nj[3] = {fullni[22].x, fullni[22].y, fullni[22].z};
        rij[0] = isqrt3; rij[1] = -isqrt3; rij[2] = isqrt3;
        E += bulk_potential(ni, nj, params, rij, 3);
    }
    if (fullni[24].pt) {
        float nj[3] = {fullni[24].x, fullni[24].y, fullni[24].z};
        rij[0] = -isqrt3; rij[1] = isqrt3; rij[2] = isqrt3;
        E += bulk_potential(ni, nj, params, rij, 3);
    }
    if (fullni[26].pt) {
        float nj[3] = {fullni[26].x, fullni[26].y, fullni[26].z};
        rij[0] = -isqrt3; rij[1] = -isqrt3; rij[2] = isqrt3;
        E += bulk_potential(ni, nj, params, rij, 3);
    }

    // --- Vertices Inferiores (Z negativo) ---
    if (fullni[21].pt) {
        float nj[3] = {fullni[21].x, fullni[21].y, fullni[21].z};
        rij[0] = isqrt3; rij[1] = isqrt3; rij[2] = -isqrt3;
        E += bulk_potential(ni, nj, params, rij, 3);
    }
    if (fullni[23].pt) {
        float nj[3] = {fullni[23].x, fullni[23].y, fullni[23].z};
        rij[0] = isqrt3; rij[1] = -isqrt3; rij[2] = -isqrt3;
        E += bulk_potential(ni, nj, params, rij, 3);
    }
    if (fullni[25].pt) {
        float nj[3] = {fullni[25].x, fullni[25].y, fullni[25].z};
        rij[0] = -isqrt3; rij[1] = isqrt3; rij[2] = -isqrt3;
        E += bulk_potential(ni, nj, params, rij, 3);
    }
    if (fullni[27].pt) {
        float nj[3] = {fullni[27].x, fullni[27].y, fullni[27].z};
        rij[0] = -isqrt3; rij[1] = -isqrt3; rij[2] = -isqrt3;
        E += bulk_potential(ni, nj, params, rij, 3);
    }
    return (float)E;
}

void GeometryStrategy::initializeBoundaries(Parameters* params) {
    // --- Configuração da Condição de Contorno em X ---
    if (strcasecmp(params->XBoundtype, "free") == 0) {
        params->XBound = &Potential::Free_Boundary;
    } else if (strcasecmp(params->XBoundtype, "periodic") == 0) {
        params->XBound = &Potential::Periodic_Boundary;
    } else {
        fprintf(stderr, "X boundary condition: %s not implemented \n", params->XBoundtype);
        exit(2);
    }

    // --- Configuração da Condição de Contorno em Y ---
    if (strcasecmp(params->YBoundtype, "free") == 0) {
        params->YBound = &Potential::Free_Boundary;
    } else if (strcasecmp(params->YBoundtype, "periodic") == 0) {
        params->YBound = &Potential::Periodic_Boundary;
    } else {
        fprintf(stderr, "Y boundary condition: %s not implemented \n", params->YBoundtype);
        exit(2);
    }

    // --- Configuração da Condição de Contorno em Z ---
    if (strcasecmp(params->ZBoundtype, "free") == 0) {
        params->ZBound = &Potential::Free_Boundary;
    } else if (strcasecmp(params->ZBoundtype, "periodic") == 0) {
        params->ZBound = &Potential::Periodic_Boundary;
    } else {
        fprintf(stderr, "Z boundary condition: %s not implemented \n", params->ZBoundtype);
        exit(2);
    }
}

// ========== Implementação da Geometria Bulk ==========

int* BulkGeometryStrategy::setPointTypes(int* pt, Parameters* params, float* surface_normals) {
    int Nx = params->Nx;
    int Ny = params->Ny;
    int Nz = params->Nz;
    
    for (int i = 0; i < Nx * Ny * Nz; i++) {
        pt[i] = 1;
    }
    return pt;
}

float BulkGeometryStrategy::calculatePotential(const nni* fullni, Parameters* params,
                                             std::vector<AnchoringStrategy*>& surfaces, 
                                             float* surface_normals) {
    float E = calculateNewmanNeighbours(fullni, params);  // 1º

    if (params->neighbourKind > 1)
        E += calculateSecondNeighbours(fullni, params);   // 2º
    if (params->neighbourKind == 3)
        E += calculateThirdNeighbours(fullni, params);    // 3º
    
    float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
    if (params->elecA != 0) 
        E += Potential::Electric_Potential(ni, params);

    return E;
}

void BulkGeometryStrategy::initializeGeometry(Parameters* params, float** surface_normals,
                                            std::vector<AnchoringStrategy*>& surfaces) {
    printf("Geometry: Bulk\n");
    
    // Alocação do vetor de normais de superfície
    *surface_normals = (float *)calloc(params->Nx * params->Ny * params->Nz * 3, sizeof(float));
    
    // Configura condições de contorno
    initializeBoundaries(params);
    
    printf("xbound  %s\n", params->XBoundtype);
    printf("ybound  %s\n", params->YBoundtype);
    printf("zbound  %s\n", params->ZBoundtype);
    printf("\n");
}

// ========== Implementação da Geometria Slab ==========

int* SlabGeometryStrategy::setPointTypes(int* pt, Parameters* params, float* surface_normals) {
    int Nx = params->Nx;
    int Ny = params->Ny;
    int Nz = params->Nz;
    
    for (int ii = 0; ii < Nx; ii++) {
        for (int jj = 0; jj < Ny; jj++) {
            int kk = 0;
            int idx_bottom = ii + Nx * (jj + Ny * kk);
            
            pt[idx_bottom] = 2; 
            surface_normals[idx_bottom * 3 + 0] = 0;
            surface_normals[idx_bottom * 3 + 1] = 0;
            surface_normals[idx_bottom * 3 + 2] = -1;

            for (kk = 1; kk < Nz - 1; kk++) {
                pt[ii + Nx * (jj + Ny * kk)] = 1; 
            }

            kk = Nz - 1;
            int idx_top = ii + Nx * (jj + Ny * kk);
            
            pt[idx_top] = 3; 
            surface_normals[idx_top * 3 + 0] = 0;
            surface_normals[idx_top * 3 + 1] = 0;
            surface_normals[idx_top * 3 + 2] = 1;
        }
    }
    return pt;
}

float SlabGeometryStrategy::calculatePotential(const nni* fullni, Parameters* params,
                                             std::vector<AnchoringStrategy*>& surfaces, 
                                             float* surface_normals) {
    float E = calculateNewmanNeighbours(fullni, params);  // 1º

    if (params->neighbourKind > 1)
        E += calculateSecondNeighbours(fullni, params);   // 2º
    if (params->neighbourKind == 3)
        E += calculateThirdNeighbours(fullni, params);    // 3º
    
    float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
    
    // Interação de Superfície - CORRIGIDO
    // ============================================================
    // MELHORIA/UNIFICAÇÃO:
    // A normal de superfície já é passada em fullni[7] (preenchida no MonteCarloStep).
    // Isso remove dependência do tipo de ponto e também corrige Sphere/Custom.
    // ============================================================
    if (fullni[0].pt > 1) {
        float s[3] = {fullni[7].x, fullni[7].y, fullni[7].z};
        E += surfaces[fullni[0].pt - 2]->calculateSurfacePotential(ni, s, params);
    }

    if (params->elecA != 0) 
        E += Potential::Electric_Potential(ni, params);

    return E;
}

void SlabGeometryStrategy::initializeGeometry(Parameters* params, float** surface_normals,
                                            std::vector<AnchoringStrategy*>& surfaces) {
    printf("Geometry: Slab\n");
    
    // Alocação do vetor de normais
    *surface_normals = (float *)calloc(params->Nx * params->Ny * params->Nz * 3, sizeof(float));
    surface_normals_ = *surface_normals;
    
    // Configura condições de contorno (Z é free no slab)
    sprintf(params->ZBoundtype, "free");
    initializeBoundaries(params);
    
    printf("xbound  %s\n", params->XBoundtype);
    printf("ybound  %s\n", params->YBoundtype);
    printf("zbound  %s\n", params->ZBoundtype);
    printf("\n");
}

// ========== Implementação da Geometria Sphere ==========

int* SphereGeometryStrategy::setPointTypes(int* pt, Parameters* params, float* surface_normals) {
    int Nx = params->Nx;
    int Ny = params->Ny;
    int Nz = params->Nz;
    float Rx, Ry, Rz;
    float Hx = Nx / 2;
    float Hy = Ny / 2;
    float Hz = Nz / 2;

    for (int ii = 0; ii < Nx; ii++) {
        for (int jj = 0; jj < Ny; jj++) {
            for (int kk = 0; kk < Nz; kk++) {
                
                Rx = ii - Hx;
                Ry = jj - Hy;
                Rz = kk - Hz;
                float radius = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);
                
                int idx = ii + Nx * (jj + Ny * kk);

                if (radius < Hz - 1) {
                    pt[idx] = 1;
                } 
                else if (radius < Hz + 1) {
                    pt[idx] = 2;
                    // OBS: radius pode ser 0 no centro, mas aqui estamos na casca (Hz±1),
                    // então não deve ocorrer. Mantemos como estava.
                    surface_normals[idx * 3 + 0] = Rx / radius;
                    surface_normals[idx * 3 + 1] = Ry / radius;
                    surface_normals[idx * 3 + 2] = Rz / radius;
                } 
                else {
                    pt[idx] = 0;
                }
            }
        }
    }
    return pt;
}

float SphereGeometryStrategy::calculatePotential(const nni* fullni, Parameters* params,
                                               std::vector<AnchoringStrategy*>& surfaces, 
                                               float* surface_normals) {
    float E = calculateNewmanNeighbours(fullni, params);  // 1º

    if (params->neighbourKind > 1)
        E += calculateSecondNeighbours(fullni, params);   // 2º
    if (params->neighbourKind == 3)
        E += calculateThirdNeighbours(fullni, params);    // 3º
    
    float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
    
    // CORREÇÃO CRÍTICA:
    // Antes estava usando normal "default temporário".
    // Agora usa a normal correta já fornecida em fullni[7].
    if (fullni[0].pt > 1) {
        float s[3] = {fullni[7].x, fullni[7].y, fullni[7].z};
        E += surfaces[fullni[0].pt - 2]->calculateSurfacePotential(ni, s, params);
    }
    
    if (params->elecA != 0) 
        E += Potential::Electric_Potential(ni, params);

    return E;
}

void SphereGeometryStrategy::initializeGeometry(Parameters* params, float** surface_normals,
                                              std::vector<AnchoringStrategy*>& surfaces) {
    printf("Geometry: Sphere\n");
    
    // Força condições de contorno livres em todas as direções
    sprintf(params->XBoundtype, "free");
    sprintf(params->YBoundtype, "free");
    sprintf(params->ZBoundtype, "free");
    
    // Alocação do vetor de normais
    *surface_normals = (float *)calloc(params->Nx * params->Ny * params->Nz * 3, sizeof(float));
    surface_normals_ = *surface_normals;
    
    initializeBoundaries(params);
    printf("\n");
}

// ========== Implementação da Geometria Custom ==========

int* CustomGeometryStrategy::setPointTypes(int* pt, Parameters* params, float* surface_normals) {
    // Implementação básica - em produção, ler de arquivo
    int Nx = params->Nx;
    int Ny = params->Ny;
    int Nz = params->Nz;
    
    for (int i = 0; i < Nx * Ny * Nz; i++) {
        pt[i] = 1; // Padrão: todos são bulk
    }
    return pt;
}

float CustomGeometryStrategy::calculatePotential(const nni* fullni, Parameters* params,
                                               std::vector<AnchoringStrategy*>& surfaces, 
                                               float* surface_normals) {
    float E = calculateNewmanNeighbours(fullni, params);  // 1º

    if (params->neighbourKind > 1)
        E += calculateSecondNeighbours(fullni, params);   // 2º
    if (params->neighbourKind == 3)
        E += calculateThirdNeighbours(fullni, params);    // 3º
    
    float ni[3] = {fullni[0].x, fullni[0].y, fullni[0].z};
    
    if (fullni[0].pt > 1) {
        float s[3] = {fullni[7].x, fullni[7].y, fullni[7].z};
        E += surfaces[fullni[0].pt - 2]->calculateSurfacePotential(ni, s, params);
    }
    
    if (params->elecA != 0) 
        E += Potential::Electric_Potential(ni, params);

    return E;
}

void CustomGeometryStrategy::initializeGeometry(Parameters* params, float** surface_normals,
                                              std::vector<AnchoringStrategy*>& surfaces) {
    printf("Geometry: Custom\n");
    
    // Alocação do vetor de normais
    *surface_normals = (float *)calloc(params->Nx * params->Ny * params->Nz * 3, sizeof(float));
    surface_normals_ = *surface_normals;
    
    initializeBoundaries(params);
    
    printf("xbound  %s\n", params->XBoundtype);
    printf("ybound  %s\n", params->YBoundtype);
    printf("zbound  %s\n", params->ZBoundtype);
    printf("\n");
}

// ========== Factory Method ==========

GeometryStrategy* GeometryStrategyFactory::create(const std::string& geometryType) {
    if (geometryType == "bulk") {
        return new BulkGeometryStrategy();
    } else if (geometryType == "slab") {
        return new SlabGeometryStrategy();
    } else if (geometryType == "sphere") {
        return new SphereGeometryStrategy();
    } else if (geometryType == "custom") {
        return new CustomGeometryStrategy();
    } else {
        std::cerr << "Unknown geometry type: " << geometryType << std::endl;
        return nullptr;
    }
}

GeometryStrategy* GeometryStrategyFactory::create(Parameters* params) {
    return create(params->geometry);
}