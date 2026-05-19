#ifndef SIMULATOR_FACTORY_H_
#define SIMULATOR_FACTORY_H_

#include "../include/parameters.h"

class Evolve;
class Geometry;

// SimulatorFactory centralizes all construction decisions that were previously
// duplicated between simulator.cpp and simulatorGPU.cpp.
//
// Usage:
//   Evolve*   ev  = SimulatorFactory::createEvolve(ni, pt, &params);
//   Geometry* geo = SimulatorFactory::createGeometry(pt, &params);
//   SimulatorFactory::setupBoundaries(params);
//   SimulatorFactory::setupPotential(params, ev);
//
// The factory reads params.use_gpu to decide between CPU (EvolveN) and
// GPU (EvolveNGPU) backends. The GPU path is only available when the
// binary is compiled with -DCUDA__ (i.e. the GPU build).
class SimulatorFactory {
public:
    // Creates the correct Evolve subclass based on params.evol and params.use_gpu.
    // Caller takes ownership of the returned pointer.
    static Evolve* createEvolve(float* ni, int* pt, Parameters* params);

    // Creates the correct Geometry subclass based on params.geometry.
    // Caller takes ownership of the returned pointer.
    static Geometry* createGeometry(int* pt, Parameters* params);

    // Resolves and assigns the boundary function pointers (XBound/YBound/ZBound).
    static void setupBoundaries(Parameters& params);

    // Resolves and assigns the bulk potential function pointer on the geometry.
    static void setupPotential(Parameters& params, Evolve* ev);
};

#endif  // SIMULATOR_FACTORY_H_
