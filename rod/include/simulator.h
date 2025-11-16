#ifndef SIMULATOR_H_
#define SIMULATOR_H_
#include <gsl/gsl_rng.h>

#include <iostream>
#include <memory> 

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
    
  std::unique_ptr<Evolve> evolve; 
  
  std::unique_ptr<int[]> pt; 
  
  ~simulator();

 private:
  std::unique_ptr<float[]> ni; 
};

#endif