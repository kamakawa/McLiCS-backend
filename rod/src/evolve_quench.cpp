#include <gsl/gsl_rng.h>
#include <math.h>
#include <omp.h>

#include <iostream>

#include "../include/define.h"
#include "../include/evolve.h"
#include "../include/geometry.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

quenchEvolveN::quenchEvolveN(float *ni, int *ppt, Parameters *params)
    // Alteração 1: A lista de inicialização Nx, Ny, Nz é removida daqui, pois é tratada na classe base EvolveN (que você corrigiu)
    : EvolveN(ni, ppt, params) { 
  printf("Initializing Step loop:\n");
  // Alteração 2: Acesso a first_file aninhado em 'mc'
  printf("Initial File Number= %d\n", params->mc.first_file);
  // Alteração 3, 4: Acesso a fn e first_file aninhados em 'mc'
  printf("Last File Number= %d\n\n", params->mc.first_file + params->mc.fn);
  this->ni = ni;
  this->pt = ppt;
  this->params = params;
};

int quenchEvolveN::run() {
  float S1, S2, sTemp;
  float tempE, E2, E;
  float vec_nt[3];
  float vec_n[3];
  float mat_n[9];
  float ang_var = 0.5;
  char fname[1000];

  int num_threads = omp_get_max_threads();
  gsl_rng **rng = (gsl_rng **)calloc(num_threads, sizeof(gsl_rng *));
  gsl_rng_env_setup();
  for (int i = 0; i < num_threads; i++) {
    rng[i] = gsl_rng_alloc(gsl_rng_ranlxs0);
    gsl_rng_set(rng[i], i);
  }

  sprintf(fname, "po.dat");
  FILE *po_file = fopen(fname, "a");
  fprintf(po_file, "ii S varS E varE\n");
  fflush(po_file);
  // Alteração 5: Acesso a T e Ti aninhados em 'potential'
  params->potential.T = params->potential.Ti;
  
  // Alteração 6, 7, 8: Acesso a MCT, MCS e fn aninhados em 'mc'
  printf("Step relaxation, for nematic molecules, using MCT=%d MCS=%d and fn=%d using %d threads\n",
         params->mc.MCT, params->mc.MCS, params->mc.fn, num_threads);
  fflush(stdout);
  
  // Alteração 9, 10, 11: Acesso a first_file e fn aninhados em 'mc'
  for (int ii = params->mc.first_file; ii < params->mc.fn + params->mc.first_file; ii++) {
    
    // Alteração 12, 13: Acesso a T e Ti aninhados em 'potential'
    params->potential.T = params->potential.Ti;
    
    // Alteração 14: Acesso a MCT aninhado em 'mc'
    for (int step = 0; step < params->mc.MCT; step++) {
      Monte_Carlo_Step(ang_var, rng);
    }
    
    // Alteração 15: Acesso a T e Tf aninhados em 'potential'
    params->potential.T = params->potential.Tf;
    
    // Alteração 16, 17: Acesso a MCT e dT aninhados em 'mc' e 'potential'
    for (int step = 0; step < params->mc.MCT*params->potential.dT; step++) {
      Monte_Carlo_Step(ang_var, rng);
    }
    S1 = 0;
    S2 = 0;
    E = 0;
    E2 = 0;
    
    // Alteração 18: Acesso a MCS aninhado em 'mc'
    for (int step = 0; step < params->mc.MCS; step++) {
      Monte_Carlo_Step(ang_var, rng);
      tempE = energy_calculator();
      E2 += tempE * tempE;
      E += tempE;
      
      // Alteração 19: Chamada de função do namespace OrderParameters
      OrderParameters::Matrice_constructor(ni, mat_n, pt, *params);
      
      // Alteração 20: Chamada de função do namespace OrderParameters
      sTemp = OrderParameters::Eigen_value_evaluation(mat_n, vec_n);
      S1 += sTemp;
      S2 += sTemp * sTemp;
    }
    
    // Alteração 21, 22: Acesso a MCS aninhado em 'mc'
    E /= params->mc.MCS;
    E2 /= params->mc.MCS;
    S1 /= params->mc.MCS;
    S2 /= params->mc.MCS;
    
    sprintf(fname, "director_field_%d.csv", ii);
    
    // Alteração 23: Chamada de função do namespace IO
    IO::print_n(fname, ni, *params, pt);
    
    fprintf(po_file, "%d %g %g %g %g\n", ii, S1, S2 - S1 * S1, E, (E2 - E * E));
    fflush(po_file);
  }

  for (int i = 0; i < num_threads; i++) gsl_rng_free(rng[i]);
  return 0;
}