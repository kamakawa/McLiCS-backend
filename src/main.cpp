#include <gsl/gsl_eigen.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <string>
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

// Removes output files from previous simulations to avoid data mixing.
static void cleanPreviousOutput() {
  printf("Cleaning up previous simulation files...\n");

  const char* fixedFiles[] = {"ic.csv", "po.dat", nullptr};
  for (int i = 0; fixedFiles[i] != nullptr; i++)
    std::filesystem::remove(fixedFiles[i]);

  for (auto& entry : std::filesystem::directory_iterator(".")) {
    const std::string name = entry.path().filename().string();
    if (name.rfind("director_field_", 0) == 0 && name.size() > 4 &&
        name.substr(name.size() - 4) == ".csv") {
      std::filesystem::remove(entry.path());
    }
  }

  printf("Cleanup done.\n\n");
}

// Reads po.dat and prints a statistical summary of the run.
static void printSimulationSummary(const Parameters* params) {
  printf("\n=== SIMULATION FINAL REPORT ===\n");

  FILE* poFile = fopen("po.dat", "r");
  if (!poFile) {
    printf("po.dat not found\n");
    printf("================================\n\n");
    return;
  }

  double sumS = 0.0, sumE = 0.0;
  double maxS = -1e9, minS = 1e9;
  double maxE = -1e9, minE = 1e9;
  int count = 0;
  char line[256];

  // Skip header
  if (!fgets(line, sizeof(line), poFile)) {
    fclose(poFile);
    printf("po.dat is empty or unreadable\n");
    return;
  }

  double T, S, varS, E, varE;
  while (fgets(line, sizeof(line), poFile)) {
    if (sscanf(line, "%lf %lf %lf %lf %lf", &T, &S, &varS, &E, &varE) == 5) {
      sumS += S;  sumE += E;
      if (S > maxS) maxS = S;
      if (S < minS) minS = S;
      if (E > maxE) maxE = E;
      if (E < minE) minE = E;
      count++;
    }
  }
  fclose(poFile);

  if (count > 0) {
    printf("Order parameter (S):\n");
    printf("  Mean:  %.4f\n", sumS / count);
    printf("  Max:   %.4f\n", maxS);
    printf("  Min:   %.4f\n", minS);
    printf("\nEnergy (E):\n");
    printf("  Mean:  %.4f\n", sumE / count);
    printf("  Max:   %.4f\n", maxE);
    printf("  Min:   %.4f\n", minE);
    printf("\nStatistics:\n");
    printf("  Saved steps: %d\n", count);
  } else {
    printf("No valid data found in po.dat\n");
  }

  printf("\nSystem configuration:\n");
  printf("  Dimensions: %dx%dx%d\n", params->Nx, params->Ny, params->Nz);
  printf("  Potential:  %s\n",        params->potential);
  printf("  Geometry:   %s\n",        params->geometry);
  printf("  Boundaries: %s / %s / %s\n",
         params->XBoundtype, params->YBoundtype, params->ZBoundtype);
  printf("================================\n\n");
}

int main(int argc, char** argv) {
  printf("### Starting McLiCS version: %s ###\n\n", MCLICS_VERSION);

  if (argc < 2) {
    printf("Error: no input file specified.\n");
    printf("Usage: %s <parameter_file>\n", argv[0]);
    return 1;
  }

  cleanPreviousOutput();

  Parameters params = read_input_file(argv[1]);
  print_parameters(params);

  char fname[1000];
  simulator* sim = new simulator(&params);
  sim->Setup_simmulation(params);

  sprintf(fname, "ic.csv");
  sim->print_n(fname, &params);
  sim->evolve->run();

  printSimulationSummary(&params);

  delete sim;
  return 0;
}
