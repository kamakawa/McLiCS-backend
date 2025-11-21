#ifndef SIMULATOR_H_
#define SIMULATOR_H_

// --- System Includes ---
#include <gsl/gsl_rng.h>
#include <iostream>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/evolve.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

class simulator {
 public:
  // Construtor e Destrutor
  simulator(Parameters *params);
  ~simulator();

  // Metodos de Configuracao e Saida
  void Setup_simmulation(Parameters &params);
  int print_n(char *fname, Parameters *params);

  // Variaveis Publicas de Estado
  Parameters *params;
  Evolve *evolve;
  int *pt;
  int Nx, Ny, Nz;

 private:
  float *ni;
};

#endif