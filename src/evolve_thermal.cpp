#include "../include/evolve.h"
#include "../include/io.h"
#include "../include/parameters.h"
#include "../include/rng_pool.h"

#include <math.h>
#include <omp.h>
#include <cstdio>

thermalEvolveN::thermalEvolveN(float* ni, int* pt, Parameters* params)
    : EvolveN(ni, pt, params) {
  printf("Initializing thermal loop:\n");
  printf("  Ti= %g\n  Tf= %g\n  dT= %g\n\n", params->Ti, params->Tf, params->dT);
}

int thermalEvolveN::run() {
  RngPool rng;
  float ang_var = 0.5f;
  char  fname[1000];

  FILE* po = open_po_file("T S varS E varE");
  const int sign = (params->dT < 0) ? -1 : 1;

  printf("Thermal sweep  Ti=%g -> Tf=%g  dT=%g  MCT=%d  MCS=%d  threads=%d\n",
         params->Ti, params->Tf, params->dT,
         params->MCT, params->MCS, rng.size());
  fflush(stdout);

  for (params->T = params->Ti;
       static_cast<int>(1e6f) * sign * (params->T - params->Tf) >= 0;
       params->T += params->dT) {

    equilibrate(params->MCT, ang_var, rng);
    auto m = measure_block(params->MCS, ang_var, rng);

    snprintf(fname, sizeof(fname),
             "director_field_%d.csv", static_cast<int>(100*(params->T + 1e-7f)));
    save_snapshot(fname);
    log_measurement(po, params->T, m);
  }

  fclose(po);
  return 0;
}
