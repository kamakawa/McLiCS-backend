#ifndef ANCHORING_STRATEGY_H_
#define ANCHORING_STRATEGY_H_

// --- System Includes ---
#include <string>
#include <vector>

// --- Project Includes ---
#include "../include/parameters.h"

// Interface Strategy para Ancoramento
class AnchoringStrategy {
public:
    virtual ~AnchoringStrategy() = default;
    
    virtual std::string getName() const = 0;
    virtual float calculateSurfacePotential(float ni[3], float s[3], Parameters* params) = 0;
    virtual void initialize(Parameters* params, int id) = 0;
    
    // Método utilitário compartilhado
    void checkParameter(bool hasValue, const std::string& parameterName, int id, const std::string& strategyName);
    
protected:
    int id_;
    Parameters* params_;
};

// Estratégias Concretas - Strong Anchoring
class StrongAnchoringStrategy : public AnchoringStrategy {
public:
    std::string getName() const override { return "Strong Anchoring"; }
    float calculateSurfacePotential(float ni[3], float s[3], Parameters* params) override;
    void initialize(Parameters* params, int id) override;
    
private:
    float W_, theta_s_, phi_s_;
};

class StrongAnchoringGHRLStrategy : public AnchoringStrategy {
public:
    std::string getName() const override { return "Strong Anchoring GHRL"; }
    float calculateSurfacePotential(float ni[3], float s[3], Parameters* params) override;
    void initialize(Parameters* params, int id) override;
    
private:
    float W_, theta_s_, phi_s_;
};

// Estratégias Concretas - Rapini-Papoular
class RPAnchoringStrategy : public AnchoringStrategy {
public:
    std::string getName() const override { return "Rapini-Papoular Anchoring"; }
    float calculateSurfacePotential(float ni[3], float s[3], Parameters* params) override;
    void initialize(Parameters* params, int id) override;
    
private:
    float W_, theta_s_, phi_s_;
};

class RPAnchoringGHRLStrategy : public AnchoringStrategy {
public:
    std::string getName() const override { return "Rapini-Papoular Anchoring GHRL"; }
    float calculateSurfacePotential(float ni[3], float s[3], Parameters* params) override;
    void initialize(Parameters* params, int id) override;
    
private:
    float W_, theta_s_, phi_s_;
};

// Estratégias Concretas - Fournier-Galatola
class FGAnchoringStrategy : public AnchoringStrategy {
public:
    std::string getName() const override { return "Fournier-Galatola like Anchoring"; }
    float calculateSurfacePotential(float ni[3], float s[3], Parameters* params) override;
    void initialize(Parameters* params, int id) override;
    
private:
    float W_;
};

class FGAnchoringGHRLStrategy : public AnchoringStrategy {
public:
    std::string getName() const override { return "Fournier-Galatola like Anchoring GHRL"; }
    float calculateSurfacePotential(float ni[3], float s[3], Parameters* params) override;
    void initialize(Parameters* params, int id) override;
    
private:
    float W_;
};

// Estratégias Concretas - Homeotropic
class HomeotropicAnchoringStrategy : public AnchoringStrategy {
public:
    std::string getName() const override { return "Homeotropic Anchoring"; }
    float calculateSurfacePotential(float ni[3], float s[3], Parameters* params) override;
    void initialize(Parameters* params, int id) override;
    
private:
    float W_;
};

class HomeotropicAnchoringGHRLStrategy : public AnchoringStrategy {
public:
    std::string getName() const override { return "Homeotropic Anchoring GHRL"; }
    float calculateSurfacePotential(float ni[3], float s[3], Parameters* params) override;
    void initialize(Parameters* params, int id) override;
    
private:
    float W_;
};

// Factory para criar estratégias de ancoramento
class AnchoringStrategyFactory {
public:
    static AnchoringStrategy* create(const std::string& anchoringType);
    static void initializeAnchoringStrategies(Parameters* params, std::vector<AnchoringStrategy*>& strategies);
};

#endif