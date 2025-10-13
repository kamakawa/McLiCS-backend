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

thermalEvolveN::thermalEvolveN(float *ni, int *ppt, Parameters *params)
    // Alteração 1: A lista de inicialização Nx, Ny, Nz é removida daqui, pois é tratada na classe base EvolveN (que você corrigiu)
    : EvolveN(ni, ppt, params) { 
  this->ni = ni;
  this->pt = ppt;
  this->params = params;
  printf("Initializing thermal loop:\n");
  // Alteração 2, 3, 4: Acesso aos parâmetros de temperatura aninhados em 'potential'
  printf("Ti= %g\n", params->potential.Ti);
  printf("Tf= %g\n", params->potential.Tf);
  printf("dT= %g\n\n", params->potential.dT);
};

int thermalEvolveN::run() {
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
  // Alteração 5: Acesso a dT aninhado em 'potential'
  int sign = -params->potential.dT / fabs(params->potential.dT);
  FILE *po_file = fopen(fname, "a");
  fprintf(po_file, "T S varS E varE\n");
  fflush(po_file);
  // Alteração 6, 7, 8: Acesso aos parâmetros de temperatura aninhados em 'potential'
  printf("Starting thermal variation, for nematic molecules, from %g to %g with step os size %g\n",
         params->potential.Ti, params->potential.Tf, params->potential.dT);
  
  // Alteração 9, 10: Acesso a MCT e MCS aninhados em 'mc'
  printf("MCT=%d MCS=%d and %d threads\n",
         params->mc.MCT, params->mc.MCS, num_threads);
  fflush(stdout);
  fflush(stdout);
  
  // Alteração 11, 12, 13, 14: Loop de variação de temperatura (T, Ti, Tf, dT)
  for (params->potential.T = params->potential.Ti; 
       (int)1e6 * sign * (params->potential.T - params->potential.Tf) >= 0; 
       params->potential.T += params->potential.dT) {
    
    // Alteração 15: Acesso a MCT aninhado em 'mc'
    for (int step = 0; step < params->mc.MCT; step++) {
      Monte_Carlo_Step(ang_var, rng);
    }
    S1 = 0;
    S2 = 0;
    E = 0;
    E2 = 0;
    // Alteração 16: Acesso a MCS aninhado em 'mc'
    for (int step = 0; step < params->mc.MCS; step++) {
      Monte_Carlo_Step(ang_var, rng);
      tempE = energy_calculator();
      E2 += tempE * tempE;
      E += tempE;
      
      // Alteração 17: Chamada de função do namespace OrderParameters
      OrderParameters::Matrice_constructor(ni, mat_n, pt, *params);
      
      // Alteração 18: Chamada de função do namespace OrderParameters
      sTemp = OrderParameters::Eigen_value_evaluation(mat_n, vec_n);
      S1 += sTemp;
      S2 += sTemp * sTemp;
    }
    // Alteração 19, 20: Acesso a MCS aninhado em 'mc'
    E /= params->mc.MCS;
    E2 /= params->mc.MCS;
    S1 /= params->mc.MCS;
    S2 /= params->mc.MCS;
    
    // Alteração 21: Acesso a T aninhado em 'potential'
    sprintf(fname, "director_field_%d.csv", (int)(100 * (params->potential.T + 1e-7)));
    
    // Alteração 22: Chamada de função do namespace IO
    IO::print_n(fname, ni, *params, pt);
    
    // Alteração 23: Acesso a T aninhado em 'potential'
    fprintf(po_file, "%g %g %g %g %g\n", params->potential.T, S1, S2 - S1 * S1, E, (E2 - E * E));
    fflush(po_file);
  }

  for (int i = 0; i < num_threads; i++)
    gsl_rng_free(rng[i]);

  return 0;
}