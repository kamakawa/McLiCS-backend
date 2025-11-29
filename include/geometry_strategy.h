#ifndef GEOMETRY_STRATEGY_H_
#define GEOMETRY_STRATEGY_H_

// --- System Includes ---
#include <vector>
#include <string>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/parameters.h"

// Forward declaration
class AnchoringStrategy;

// Interface Strategy para Geometria
class GeometryStrategy {
public:
    virtual ~GeometryStrategy() = default;
    
    // Métodos principais
    virtual int* setPointTypes(int* pt, Parameters* params, float* surface_normals) = 0;
    virtual float calculatePotential(const nni fullni[7], Parameters* params, 
                                   std::vector<AnchoringStrategy*>& surfaces, 
                                   float* surface_normals) = 0;
    virtual void initializeGeometry(Parameters* params, float** surface_normals, 
                                  std::vector<AnchoringStrategy*>& surfaces) = 0;
    
    // Getters
    virtual std::string getName() const = 0;
    virtual int getNumSurfaces() const = 0;
    virtual float* getSurfaceNormals() = 0;
    
    // NOVO: Setter para potencial bulk
    virtual void setBulkPotential(float (*potential)(float ni[3], float nj[3], Parameters* params, float rij[3], int nk)) {
        bulk_potential = potential;
    }
    
    // Métodos utilitários compartilhados
    float calculateNewmanNeighbours(const nni fullni[], Parameters* params);
    float calculateSecondNeighbours(const nni fullni[], Parameters* params);
    float calculateThirdNeighbours(const nni fullni[], Parameters* params);
    void initializeBoundaries(Parameters* params);
    
protected:
    // Ponteiro para função de potencial bulk
    float (*bulk_potential)(float ni[3], float nj[3], Parameters* params, float rij[3], int nk) = nullptr;
};

// Estratégias Concretas
class BulkGeometryStrategy : public GeometryStrategy {
public:
    int* setPointTypes(int* pt, Parameters* params, float* surface_normals) override;
    float calculatePotential(const nni fullni[7], Parameters* params,
                           std::vector<AnchoringStrategy*>& surfaces, 
                           float* surface_normals) override;
    void initializeGeometry(Parameters* params, float** surface_normals,
                          std::vector<AnchoringStrategy*>& surfaces) override;
    
    std::string getName() const override { return "Bulk Geometry"; }
    int getNumSurfaces() const override { return 0; }
    float* getSurfaceNormals() override { return nullptr; }
};

class SlabGeometryStrategy : public GeometryStrategy {
public:
    int* setPointTypes(int* pt, Parameters* params, float* surface_normals) override;
    float calculatePotential(const nni fullni[7], Parameters* params,
                           std::vector<AnchoringStrategy*>& surfaces, 
                           float* surface_normals) override;
    void initializeGeometry(Parameters* params, float** surface_normals,
                          std::vector<AnchoringStrategy*>& surfaces) override;
    
    std::string getName() const override { return "Slab Geometry"; }
    int getNumSurfaces() const override { return 2; }
    float* getSurfaceNormals() override { return surface_normals_; }
    
private:
    float* surface_normals_ = nullptr;
};

class SphereGeometryStrategy : public GeometryStrategy {
public:
    int* setPointTypes(int* pt, Parameters* params, float* surface_normals) override;
    float calculatePotential(const nni fullni[7], Parameters* params,
                           std::vector<AnchoringStrategy*>& surfaces, 
                           float* surface_normals) override;
    void initializeGeometry(Parameters* params, float** surface_normals,
                          std::vector<AnchoringStrategy*>& surfaces) override;
    
    std::string getName() const override { return "Sphere Geometry"; }
    int getNumSurfaces() const override { return 1; }
    float* getSurfaceNormals() override { return surface_normals_; }
    
private:
    float* surface_normals_ = nullptr;
};

class CustomGeometryStrategy : public GeometryStrategy {
public:
    int* setPointTypes(int* pt, Parameters* params, float* surface_normals) override;
    float calculatePotential(const nni fullni[7], Parameters* params,
                           std::vector<AnchoringStrategy*>& surfaces, 
                           float* surface_normals) override;
    void initializeGeometry(Parameters* params, float** surface_normals,
                          std::vector<AnchoringStrategy*>& surfaces) override;
    
    std::string getName() const override { return "Custom Geometry"; }
    int getNumSurfaces() const override { return nSurfaces_; }
    float* getSurfaceNormals() override { return surface_normals_; }
    
private:
    float* surface_normals_ = nullptr;
    int nSurfaces_ = 0;
};

// Factory para criar estratégias de geometria
class GeometryStrategyFactory {
public:
    static GeometryStrategy* create(const std::string& geometryType);
    static GeometryStrategy* create(Parameters* params);
};

#endif