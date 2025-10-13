#include <gsl/gsl_eigen.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>

#include <iostream>

#include "../include/define.h"
#include "../include/evolve.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"
#include "../include/simulator.h"

#define MCLICS_VERSION "0.1"
int main(int argc, char **argv) {
  printf("### Starting McLiCS version: %s ###\n\n",MCLICS_VERSION);
  
  // Alteração 1: Chamada de função do namespace IO
  Parameters params = IO::read_input_file(argv[1]);
  
  // Alteração 2: Chamada de função do namespace IO
  IO::print_parameters(params);
  
  char fname[1000];
  simulator *sim = new simulator(&params);
  sim->Setup_simmulation(params);

  sprintf(fname, "ic.csv");
  sim->print_n(fname, &params);
  sim->evolve->run();

  return 0;
}