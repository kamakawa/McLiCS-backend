#include "../include/evolve.h"
#include "../include/io.h"
#include "../include/parameters.h"
#include "../include/rng_pool.h"

#include <math.h>
#include <omp.h>
#include <cstdio>

stepEvolveN::stepEvolveN(float* ni, int* pt, Parameters* params)
    : EvolveN(ni, pt, params) {
  printf("Initializing step loop:\n");
  printf("  First file= %d  Last file= %d\n\n",
         params->first_file, params->first_file + params->fn);
}

int stepEvolveN::run() {
  RngPool rng;
  float ang_var = 0.5f;
  char  fname[1000];

  FILE* po = open_po_file("ii S varS E varE");
  params->T = params->Ti;

  printf("Step relaxation  MCT=%d  MCS=%d  fn=%d  threads=%d\n",
         params->MCT, params->MCS, params->fn, rng.size());
  fflush(stdout);

  for (int ii = params->first_file; ii < params->fn + params->first_file; ii++) {
    equilibrate(params->MCT, ang_var, rng);
    auto m = measure_block(params->MCS, ang_var, rng);

    snprintf(fname, sizeof(fname), "director_field_%d.csv", ii);
    save_snapshot(fname);
    log_measurement(po, static_cast<float>(ii), m);
  }

  fclose(po);
  return 0;
}
