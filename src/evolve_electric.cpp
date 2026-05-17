#include "../include/evolve.h"
#include "../include/io.h"
#include "../include/parameters.h"
#include "../include/rng_pool.h"

#include <math.h>
#include <omp.h>
#include <cstdio>

electricEvolveN::electricEvolveN(float* ni, int* pt, Parameters* params)
    : EvolveN(ni, pt, params) {
  printf("Initializing electric loop:\n");
  printf("  Ei= %g\n  Ef= %g\n  dE= %g\n\n",
         params->elecEi, params->elecEf, params->elecdE);
}

int electricEvolveN::run() {
  RngPool rng;
  float ang_var = 0.5f;
  char  fname[1000];

  FILE* po = open_po_file("E S varS Epot varEpot");
  params->T = params->Ti;
  const int sign = (params->elecdE < 0) ? -1 : 1;

  printf("Electric field sweep  Ei=%g -> Ef=%g  dE=%g  MCT=%d  MCS=%d  threads=%d\n",
         params->elecEi, params->elecEf, params->elecdE,
         params->MCT, params->MCS, rng.size());
  fflush(stdout);

  for (params->elecE = params->elecEi;
       static_cast<int>(1e6f) * sign * (params->elecE - params->elecEf) >= 0;
       params->elecE += params->elecdE) {

    equilibrate(params->MCT, ang_var, rng);
    auto m = measure_block(params->MCS, ang_var, rng);

    snprintf(fname, sizeof(fname),
             "director_field_%d.csv", static_cast<int>(100*(params->elecE + 1e-7f)));
    save_snapshot(fname);
    log_measurement(po, params->elecE, m);
  }

  fclose(po);
  return 0;
}
