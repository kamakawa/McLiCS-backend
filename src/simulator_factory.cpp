#include "../include/simulator_factory.h"

#include <cstdio>
#include <cstdlib>
#include <strings.h>

#include "../include/evolve.h"
#include "../include/geometry.h"
#include "../include/potential.h"

#ifdef CUDA__
#include "../include/evolve.cuh"
#endif

// ---------------------------------------------------------------------------
// Evolve factory
// ---------------------------------------------------------------------------
Evolve* SimulatorFactory::createEvolve(float* ni, int* pt, Parameters* params) {
    const char* evol = params->evol;

#ifdef CUDA__
    if (params->use_gpu) {
        if      (strcasecmp(evol, "thermal") == 0) { std::printf("  Evolution    : thermal (GPU)\n");  return new thermalEvolveNGPU(ni, pt, params); }
        else if (strcasecmp(evol, "step")    == 0) { std::printf("  Evolution    : step (GPU)\n");     return new stepEvolveNGPU(ni, pt, params); }
        else if (strcasecmp(evol, "quench")  == 0) { std::printf("  Evolution    : quench (GPU)\n");   return new quenchEvolveNGPU(ni, pt, params); }
        else if (strcasecmp(evol, "electric")== 0) { std::printf("  Evolution    : electric (GPU)\n"); return new electricEvolveNGPU(ni, pt, params); }
        else {
            std::fprintf(stderr, "Error: evol '%s' not recognised. Options: thermal | step | quench | electric\n", evol);
            std::exit(2);
        }
    }
#else
    if (params->use_gpu) {
        std::fprintf(stderr,
            "Error: device=gpu requested but this binary was compiled without CUDA.\n"
            "       Rebuild with 'make' (GPU target) or change 'device cpu' in the parameter file.\n");
        std::exit(2);
    }
#endif

    // CPU path
    if      (strcasecmp(evol, "thermal") == 0) { std::printf("  Evolution    : thermal (CPU)\n");  return new thermalEvolveN(ni, pt, params); }
    else if (strcasecmp(evol, "step")    == 0) { std::printf("  Evolution    : step (CPU)\n");     return new stepEvolveN(ni, pt, params); }
    else if (strcasecmp(evol, "quench")  == 0) { std::printf("  Evolution    : quench (CPU)\n");   return new quenchEvolveN(ni, pt, params); }
    else if (strcasecmp(evol, "electric")== 0) { std::printf("  Evolution    : electric (CPU)\n"); return new electricEvolveN(ni, pt, params); }
    else {
        std::fprintf(stderr, "Error: evol '%s' not recognised. Options: thermal | step | quench | electric\n", evol);
        std::exit(2);
    }
}

// ---------------------------------------------------------------------------
// Geometry factory
// ---------------------------------------------------------------------------
Geometry* SimulatorFactory::createGeometry(int* pt, Parameters* params) {
    const char* geo = params->geometry;

    if      (strcasecmp(geo, "bulk")   == 0) return new Bulk_Geometry(pt, params);
    else if (strcasecmp(geo, "slab")   == 0) return new Slab_Geometry(pt, params);
    else if (strcasecmp(geo, "sphere") == 0) return new Sphere_Geometry(pt, params);
    else if (strcasecmp(geo, "custom") == 0) return new Custom_Geometry(pt, params);
    else {
        std::fprintf(stderr, "Error: geometry '%s' not recognised. Options: bulk | slab | sphere | custom\n", geo);
        std::exit(2);
    }
}

// ---------------------------------------------------------------------------
// Boundary conditions
// ---------------------------------------------------------------------------
void SimulatorFactory::setupBoundaries(Parameters& p) {
    auto resolve = [](const char* name, const char* axis) -> int(*)(int&, int) {
        if      (strcasecmp(name, "free")     == 0) return &Free_Boundary;
        else if (strcasecmp(name, "periodic") == 0) return &Periodic_Boundary;
        else {
            std::fprintf(stderr, "Error: %s boundary '%s' not implemented. Options: free | periodic\n", axis, name);
            std::exit(2);
        }
    };

    p.XBound = resolve(p.XBoundtype, "X");
    p.YBound = resolve(p.YBoundtype, "Y");
    p.ZBound = resolve(p.ZBoundtype, "Z");
}

// ---------------------------------------------------------------------------
// Bulk potential
// ---------------------------------------------------------------------------
void SimulatorFactory::setupPotential(Parameters& p, Evolve* ev) {
    const char* pot = p.potential;

    if (strcasecmp(pot, "ll") * strcasecmp(pot, "lebwohl-lasher") == 0) {
        ev->geometry->bulk_potential = &Bulk_Energy_Lebwohl_Lasher;
        std::printf("  Potential    : Lebwohl-Lasher\n");
    } else if (strcasecmp(pot, "ghrl") * strcasecmp(pot, "grun-hess") == 0) {
        ev->geometry->bulk_potential = &Bulk_Energy_GHRL;
        setGHRL(p);
        std::printf("  Potential    : Gruhn-Hess (GHRL)\n");
    } else if (strcasecmp(pot, "pear") == 0) {
        ev->geometry->bulk_potential = &Bulk_Energy_Selinger_Pear;
        std::printf("  Potential    : Splay-bend (Pear)\n");
    } else {
        std::fprintf(stderr, "Error: potential '%s' not implemented. Options: ll | ghrl | pear\n", pot);
        std::exit(2);
    }
}
