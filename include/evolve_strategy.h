#ifndef EVOLVE_STRATEGY_H_
#define EVOLVE_STRATEGY_H_

#include <gsl/gsl_rng.h>
#include <vector>
#include "../include/define.h"
#include "../include/parameters.h"

// Forward declarations
class GeometryStrategy;
class AnchoringStrategy;

// Interface Strategy para Evolução
class EvolveStrategy {
public:
    virtual ~EvolveStrategy() = default;
    
    virtual int run(float* ni, int* pt, Parameters* params, 
                   GeometryStrategy* geometry_strategy,
                   float* surface_normals,
                   std::vector<AnchoringStrategy*>& anchoring_strategies) = 0;
    
    virtual std::string getName() const = 0;
    
protected:
    // Métodos utilitários compartilhados - ATUALIZADOS
    void monteCarloStep(float& ang_var, gsl_rng** r, float* ni, int* pt, 
                       Parameters* params, GeometryStrategy* geometry_strategy,
                       float* surface_normals,
                       std::vector<AnchoringStrategy*>& anchoring_strategies);
    
    float energyCalculator(float* ni, int* pt, Parameters* params,
                          GeometryStrategy* geometry_strategy,
                          float* surface_normals,
                          std::vector<AnchoringStrategy*>& anchoring_strategies);
    
    void initializeRNG(gsl_rng*** rng, int num_threads);
    void cleanupRNG(gsl_rng** rng, int num_threads);
};

// Estratégia de Evolução Térmica
class ThermalEvolveStrategy : public EvolveStrategy {
public:
    int run(float* ni, int* pt, Parameters* params,
           GeometryStrategy* geometry_strategy,
           float* surface_normals,
           std::vector<AnchoringStrategy*>& anchoring_strategies) override;
    
    std::string getName() const override { return "Thermal Evolution"; }
    
private:
    void thermalLoop(float* ni, int* pt, Parameters* params,
                    GeometryStrategy* geometry_strategy,
                    float* surface_normals,
                    std::vector<AnchoringStrategy*>& anchoring_strategies,
                    gsl_rng** rng, FILE* po_file);
};

// Estratégia de Evolução por Passos
class StepEvolveStrategy : public EvolveStrategy {
public:
    int run(float* ni, int* pt, Parameters* params,
           GeometryStrategy* geometry_strategy,
           float* surface_normals,
           std::vector<AnchoringStrategy*>& anchoring_strategies) override;
    
    std::string getName() const override { return "Step Evolution"; }
};

// Estratégia de Quench (Resfriamento Rápido)
class QuenchEvolveStrategy : public EvolveStrategy {
public:
    int run(float* ni, int* pt, Parameters* params,
           GeometryStrategy* geometry_strategy,
           float* surface_normals,
           std::vector<AnchoringStrategy*>& anchoring_strategies) override;
    
    std::string getName() const override { return "Quench Evolution"; }
    
private:
    void quenchLoop(float* ni, int* pt, Parameters* params,
                   GeometryStrategy* geometry_strategy,
                   float* surface_normals,
                   std::vector<AnchoringStrategy*>& anchoring_strategies,
                   gsl_rng** rng, FILE* po_file);
};

// Estratégia de Evolução com Campo Elétrico
class ElectricEvolveStrategy : public EvolveStrategy {
public:
    int run(float* ni, int* pt, Parameters* params,
           GeometryStrategy* geometry_strategy,
           float* surface_normals,
           std::vector<AnchoringStrategy*>& anchoring_strategies) override;
    
    std::string getName() const override { return "Electric Field Evolution"; }
    
private:
    void electricLoop(float* ni, int* pt, Parameters* params,
                     GeometryStrategy* geometry_strategy,
                     float* surface_normals,
                     std::vector<AnchoringStrategy*>& anchoring_strategies,
                     gsl_rng** rng, FILE* po_file);
};

// Factory para criar estratégias
class EvolveStrategyFactory {
public:
    static EvolveStrategy* create(const std::string& evolveType);
    static EvolveStrategy* create(Parameters* params);
};

#endif