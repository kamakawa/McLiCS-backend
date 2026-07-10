#include <dirent.h>
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

// ─── Helpers ────────────────────────────────────────────────────────────────

static void print_separator(const char *seg = "=", int count = 60) {
  for (int i = 0; i < count; i++) fputs(seg, stdout);
  putchar('\n');
}

static void print_banner() {
  print_separator("=");
  printf("  McLiCS - Monte Carlo Liquid Crystal Simulator (v%s)\n", MCLICS_VERSION);
  print_separator("=");
  printf("\n");
}

// ─── cleanPreviousOutput ────────────────────────────────────────────────────
// Removes leftover output files (ic.csv, po.dat, director_field_*.csv) from a
// previous run so results never get appended to old data.

static bool has_suffix(const char *name, const char *suffix) {
  size_t nameLen = strlen(name);
  size_t sufLen = strlen(suffix);
  if (sufLen > nameLen) return false;
  return strcmp(name + (nameLen - sufLen), suffix) == 0;
}

static bool has_prefix(const char *name, const char *prefix) {
  return strncmp(name, prefix, strlen(prefix)) == 0;
}

static void cleanPreviousOutput() {
  printf("  Removing output files from previous run... ");
  fflush(stdout);
  int removed = 0;
  DIR *dir = opendir(".");
  if (dir) {
    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
      const char *name = entry->d_name;
      if (strcmp(name, "ic.csv") == 0 || strcmp(name, "po.dat") == 0 ||
          (has_prefix(name, "director_field_") && has_suffix(name, ".csv"))) {
        remove(name);
        removed++;
      }
    }
    closedir(dir);
  }
  printf("%d file(s) removed.\n\n", removed);
}

int main(int argc, char **argv) {
  print_banner();

  cleanPreviousOutput();

  Parameters params = read_input_file(argc > 1 ? argv[1] : NULL);
  print_parameters(params);

  char fname[1000];
  simulator *sim = new simulator(&params);
  sim->Setup_simmulation(params);

  sprintf(fname, "ic.csv");
  sim->print_n(fname, &params);
  sim->evolve->run();

  print_separator("=");
  printf("  SIMULATION COMPLETE\n");
  print_separator("=");
  printf("\n");

  return 0;
}
