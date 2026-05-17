#ifndef EVOL_H_
#define EVOL_H_
#include <iostream>
#include <cstdio>

#include "../include/define.h"
#include "../include/geometry.h"
#include "../include/parameters.h"
#include "../include/potential.h"
#include "../include/rng_pool.h"

#ifdef CUDA__
#include <cuda_runtime.h>
#include <curand_kernel.h>
#else
struct dim3 {
  unsigned int x;
  unsigned int y;
  unsigned int z;
};
#endif

// ---------------------------------------------------------------------------
// Measurement — carries the averaged observables from one MCS block.
// ---------------------------------------------------------------------------
struct Measurement {
  float S1  = 0.f;  // <S>
  float S2  = 0.f;  // <S^2>
  float E   = 0.f;  // <E>
  float E2  = 0.f;  // <E^2>
};

// ---------------------------------------------------------------------------
// Base class
// ---------------------------------------------------------------------------
class Evolve {
 public:
  virtual ~Evolve() = default;

  virtual int   run()             { return 0; }
  virtual float latice_Potential(){ return 0; }

  unsigned int Nx = 0, Ny = 0, Nz = 0;
  Geometry*    geometry    = nullptr;
  int          VallidPoints = 0;

  void check_Points(int* pt, Parameters params);
};

// ---------------------------------------------------------------------------
// EvolveN — CPU Monte Carlo base with shared helpers
// ---------------------------------------------------------------------------
class EvolveN : public Evolve {
 public:
  EvolveN(float* ni, int* pt, Parameters* params);
  ~EvolveN() override;

  virtual int run() { return 0; }

  // Core MC primitives
  void  Monte_Carlo_Step(float& ang_var, gsl_rng** r);
  float energy_calculator();

  // Template-Method helpers — called by subclass run() implementations
  // to avoid duplicating the measurement loop.
  void        equilibrate(int mct, float& ang_var, RngPool& rng);
  Measurement measure_block(int mcs, float& ang_var, RngPool& rng);

  // Utility: open po.dat and write header if file is empty.
  static FILE* open_po_file(const char* header);
  // Utility: save snapshot and append one data row to po.dat.
  void save_snapshot(const char* fname);
  void log_measurement(FILE* po, float key, const Measurement& m);

 protected:
  dim3 tick{};
#ifdef CUDA__
  curandState*  d_rngStates = nullptr;
#endif
  int*          d_pt     = nullptr;
  float*        d_T      = nullptr;
  unsigned int* d_acc    = nullptr;
  Parameters*   d_params = nullptr;

  int    Nx, Ny, Nz;
  float* ni;
  int*   pt;
  Parameters* params;

  // GSL eigen workspace — allocated once, reused every call (thread-safe
  // because measure_block is called from the main thread only).
  gsl_eigen_symmv_workspace* eigen_ws_  = nullptr;
  gsl_vector*                eigen_eval_ = nullptr;
  gsl_matrix*                eigen_evec_ = nullptr;
  gsl_matrix*                eigen_m_    = nullptr;
};

// ---------------------------------------------------------------------------
// CPU subclasses — only override run()
// ---------------------------------------------------------------------------
class thermalEvolveN : public EvolveN {
 public:
  thermalEvolveN(float* ni, int* pt, Parameters* params);
  int run() override;
};

class stepEvolveN : public EvolveN {
 public:
  stepEvolveN(float* ni, int* pt, Parameters* params);
  int run() override;
};

class quenchEvolveN : public EvolveN {
 public:
  quenchEvolveN(float* ni, int* pt, Parameters* params);
  int run() override;
};

class electricEvolveN : public EvolveN {
 public:
  electricEvolveN(float* ni, int* pt, Parameters* params);
  int run() override;
};

#endif
