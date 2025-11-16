#include <gsl/gsl_rng.h>

#include <iostream>

#include "../include/define.h"
#include "../include/evolve.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

class simulator {
 public:
  simulator(Parameters *params);
  Parameters *params;
  void Setup_simmulation(Parameters &params);
  int print_n(char *fname, Parameters *params);
  int Nx, Ny, Nz;
  //~ virtual void Setup_simmulation(parameters &params){};
  //~ virtual int Evol(gsl_rng * rng, parameters *params){return 0;};
  Evolve *evolve;
  int *pt;
  ~simulator();

 private:
  float *ni;
};
