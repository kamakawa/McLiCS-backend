#include "../include/simulator.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "../include/evolve.h"
#include "../include/geometry.h"
#include "../include/ic.h"
#include "../include/io.h"
#include "../include/parameter_order.h"
#include "../include/potential.h"
#include "../include/simulator_factory.h"

simulator::simulator(Parameters* params_in)
    : params(params_in),
      Nx(params_in ? params_in->Nx : 0),
      Ny(params_in ? params_in->Ny : 0),
      Nz(params_in ? params_in->Nz : 0) {}

simulator::~simulator() {
    delete evolve;
    evolve = nullptr;

    std::free(ni);
    ni = nullptr;

    std::free(pt);
    pt = nullptr;
}

void simulator::Setup_simmulation(Parameters& p) {
    delete evolve;
    evolve = nullptr;
    std::free(ni); ni = nullptr;
    std::free(pt); pt = nullptr;

    const int Np = p.Nx * p.Ny * p.Nz;
    const int Nn = 3 * Np;

    ni = static_cast<float*>(std::calloc(static_cast<std::size_t>(Nn), sizeof(float)));
    pt = static_cast<int*>  (std::calloc(static_cast<std::size_t>(Np), sizeof(int)));

    if (!ni || !pt) {
        std::fprintf(stderr, "Out of memory allocating director/point-type arrays "
                             "(Nx=%d Ny=%d Nz=%d)\n", p.Nx, p.Ny, p.Nz);
        std::exit(2);
    }

    evolve           = SimulatorFactory::createEvolve(ni, pt, &p);
    evolve->geometry = SimulatorFactory::createGeometry(pt, &p);

    SimulatorFactory::setupBoundaries(p);
    evolve->geometry->Boundary_Init(&p);
    evolve->check_Points(pt, p);
    apply_Initial_Condidions(ni, pt, p);
    SimulatorFactory::setupPotential(p, evolve);
}

int simulator::print_n(const char* fname, const Parameters* p) const {
    if (!fname || !p || !ni || !pt) {
        std::fprintf(stderr, "print_n: called with null data\n");
        return 1;
    }

    FILE* output = std::fopen(fname, "w");
    if (!output) {
        std::perror(fname);
        return 1;
    }

    std::fprintf(output, "x,y,z,nx,ny,nz,S,pt\n");
    for (int k = 0; k < p->Nz; k++) {
        for (int j = 0; j < p->Ny; j++) {
            for (int i = 0; i < p->Nx; i++) {
                float S = pti(i, j, k)
                    ? lattice_order_parameter(ni, pt, i, j, k, *p)
                    : 1.0f;
                std::fprintf(output, "%d,%d,%d,%g,%g,%g,%g,%d\n",
                             i, j, k,
                             nix(i, j, k), niy(i, j, k), niz(i, j, k),
                             S,
                             pti(i, j, k));
            }
        }
    }

    std::fclose(output);
    std::printf("  Snapshot -> %s\n", fname);
    return 0;
}
