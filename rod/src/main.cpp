#include <gsl/gsl_eigen.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>

#include <iostream>
#include <memory> // Necessário para std::unique_ptr
#include <string> // Necessário para std::string

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
  
  if (argc < 2) {
      fprintf(stderr, "Uso: %s <input_file>\n", argv[0]);
      return 1;
  }

  Parameters params = IO::read_input_file(argv[1]);
  
  IO::print_parameters(params);
  
  std::unique_ptr<simulator> sim = std::make_unique<simulator>(&params);
  
  sim->Setup_simmulation(params);

  const std::string fname = "ic.csv"; 
  
  sim->print_n(fname, params); 
  
  sim->evolve->run();

  return 0;
}