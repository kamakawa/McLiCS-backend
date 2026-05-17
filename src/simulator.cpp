#include <cstdlib>
#include <cstring>
#include <strings.h>
#include <iostream>

#include "../include/ic.h"
#include "../include/simulator.h"
#include "../include/evolve.h"
#include "../include/geometry.h"
#include "../include/potential.h"
#include "../include/io.h"

// GPU-specific headers (only available when CUDA__ is defined)
#ifdef CUDA__
#include "../include/evolve.cuh"
#endif

// ---------------------------------------------------------------------------
// Internal factory helpers — isolate all if-else chains here so that
// Setup_simmulation stays readable.
// ---------------------------------------------------------------------------

static Evolve* makeEvolve(float* ni, int* pt, Parameters* p) {
#ifdef CUDA__
  if (strcasecmp(p->evol, "thermal") == 0 || strcasecmp(p->evol, "thermalGPU") == 0) {
    std::printf("[CUDA] Using thermal evolution on GPU\n");
    return new thermalEvolveNGPU(ni, pt, p);
  }
  if (strcasecmp(p->evol, "step") == 0 || strcasecmp(p->evol, "stepGPU") == 0) {
    std::printf("[CUDA] Using step evolution on GPU\n");
    return new stepEvolveNGPU(ni, pt, p);
  }
  if (strcasecmp(p->evol, "quench") == 0 || strcasecmp(p->evol, "quenchGPU") == 0) {
    std::printf("[CUDA] Using quench evolution on GPU\n");
    return new quenchEvolveNGPU(ni, pt, p);
  }
  if (strcasecmp(p->evol, "electric") == 0 || strcasecmp(p->evol, "electricGPU") == 0) {
    std::printf("[CUDA] Using electric evolution on GPU\n");
    return new electricEvolveNGPU(ni, pt, p);
  }
#else
  if (strcasecmp(p->evol, "thermal") == 0) return new thermalEvolveN(ni, pt, p);
  if (strcasecmp(p->evol, "step")    == 0) return new stepEvolveN(ni, pt, p);
  if (strcasecmp(p->evol, "quench")  == 0) return new quenchEvolveN(ni, pt, p);
  if (strcasecmp(p->evol, "electric")== 0) return new electricEvolveN(ni, pt, p);
#endif
  std::fprintf(stderr,
    "Evolve '%s' not defined. Try: thermal | step | quench | electric\n", p->evol);
  std::exit(2);
}

static Geometry* makeGeometry(int* pt, Parameters* p) {
  if (strcasecmp(p->geometry, "bulk")   == 0) return new Bulk_Geometry(pt, p);
  if (strcasecmp(p->geometry, "slab")   == 0) return new Slab_Geometry(pt, p);
  if (strcasecmp(p->geometry, "sphere") == 0) return new Sphere_Geometry(pt, p);
  if (strcasecmp(p->geometry, "custom") == 0) return new Custom_Geometry(pt, p);
  std::fprintf(stderr,
    "Geometry '%s' not defined. Try: bulk | slab | sphere | custom\n", p->geometry);
  std::exit(2);
}

static void applyBoundary(Parameters& p) {
  auto pick = [](const char* name, const char* axis) -> int(*)(int&,int) {
    if (strcasecmp(name, "free")     == 0) return &Free_Boundary;
    if (strcasecmp(name, "periodic") == 0) return &Periodic_Boundary;
    std::fprintf(stderr, "%s boundary condition '%s' not implemented\n", axis, name);
    std::exit(2);
  };
  p.XBound = pick(p.XBoundtype, "X");
  p.YBound = pick(p.YBoundtype, "Y");
  p.ZBound = pick(p.ZBoundtype, "Z");
}

static void applyBulkPotential(Evolve* ev, Parameters& p) {
  if (strcasecmp(p.potential, "ll") * strcasecmp(p.potential, "lebwohl-lahser") == 0) {
    ev->geometry->bulk_potential = &Bulk_Energy_Lebwohl_Lasher;
    std::printf("Using Lebwohl-Lasher potential\n");
  } else if (strcasecmp(p.potential, "ghrl") * strcasecmp(p.potential, "grun-hess") == 0) {
    ev->geometry->bulk_potential = &Bulk_Energy_GHRL;
    setGHRL(p);
    std::printf("Using Gruhn-Hess (GHRL) potential\n");
  } else if (strcasecmp(p.potential, "pear") == 0) {
    ev->geometry->bulk_potential = &Bulk_Energy_Selinger_Pear;
    std::printf("Using splay-bend (Pear) potential\n");
  } else {
    std::fprintf(stderr,
      "Potential '%s' not implemented. Try: LL | GHRL | pear\n", p.potential);
    std::exit(2);
  }
}

// ---------------------------------------------------------------------------
// simulator implementation
// ---------------------------------------------------------------------------

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
  // Release any previous state
  delete evolve;
  evolve = nullptr;
  std::free(ni); ni = nullptr;
  std::free(pt); pt = nullptr;

  const int Np = p.Nx * p.Ny * p.Nz;
  const int Nn = 3 * Np;

  ni = static_cast<float*>(std::calloc(static_cast<size_t>(Nn), sizeof(float)));
  pt = static_cast<int*>  (std::calloc(static_cast<size_t>(Np), sizeof(int)));

  if (!ni || !pt) {
    std::fprintf(stderr,
      "Out of memory allocating ni/pt (Nx=%d Ny=%d Nz=%d)\n", p.Nx, p.Ny, p.Nz);
    std::exit(2);
  }

  evolve           = makeEvolve(ni, pt, &p);
  evolve->geometry = makeGeometry(pt, &p);

  applyBoundary(p);

  evolve->geometry->Boundary_Init(&p);
  evolve->check_Points(pt, p);
  apply_Initial_Condidions(ni, pt, p);

  applyBulkPotential(evolve, p);
}

int simulator::print_n(const char* fname, const Parameters* p) const {
  if (!fname || !p || !ni || !pt) {
    std::fprintf(stderr, "print_n called with null data\n");
    return 1;
  }

  FILE* output = std::fopen(fname, "w");
  if (!output) {
    std::perror(fname);
    return 1;
  }

  std::fprintf(output, "x,y,z,nx,ny,nz,s,pt\n");
  for (int k = 0; k < p->Nz; k++) {
    for (int j = 0; j < p->Ny; j++) {
      for (int i = 0; i < p->Nx; i++) {
        std::fprintf(output, "%d,%d,%d,%g,%g,%g,1,%d\n",
                     i, j, k,
                     nix(i, j, k), niy(i, j, k), niz(i, j, k),
                     pti(i, j, k));
      }
    }
  }

  std::fclose(output);
  std::printf("Snapshot saved in %s\n", fname);
  return 0;
}
