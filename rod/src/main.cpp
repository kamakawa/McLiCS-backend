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
  printf("### Starting McLiCS version: %s ###\n\n", MCLICS_VERSION);
  
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
    return 1;
  }
  
  // Leitura dos parâmetros de entrada
  Parameters params = read_input_file(argv[1]);
  print_parameters(params);
  
  // Preparação do nome do arquivo de saída inicial
  char fname[1000];
  sprintf(fname, "ic.csv");
  
  // Criação e configuração do simulador
  simulator *sim = new simulator(&params);
  sim->Setup_simmulation(params);
  
  // Salva condição inicial e executa simulação
  sim->print_n(fname, &params);
  sim->evolve->run();
  
  // CORREÇÃO: Liberar memória do simulador
  delete sim;
  
  return 0;
}