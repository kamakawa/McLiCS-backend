#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <filesystem>
#include <string>

#include "../include/define.h"
#include "../include/evolve.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"
#include "../include/simulator.h"

#define MCLICS_VERSION "0.2"

namespace fs = std::filesystem;

// ─── Helpers ────────────────────────────────────────────────────────────────

static void print_separator(const char *seg = "─", int count = 60) {
  for (int i = 0; i < count; i++) fputs(seg, stdout);
  putchar('\n');
}

static void print_banner() {
  print_separator("═");
  printf("  McLiCS v%s  —  Monte Carlo Liquid Crystal Simulator\n", MCLICS_VERSION);
  print_separator("═");
  printf("\n");
}

// ─── cleanPreviousOutput ────────────────────────────────────────────────────

static void cleanPreviousOutput() {
  printf("  Removing output files from previous run... ");
  fflush(stdout);
  int removed = 0;
  for (const auto& entry : fs::directory_iterator(".")) {
    const std::string name = entry.path().filename().string();
    if (name == "ic.csv" || name == "po.dat" ||
        (name.rfind("director_field_", 0) == 0 && name.size() >= 4 &&
         name.substr(name.size() - 4) == ".csv")) {
      fs::remove(entry.path());
      removed++;
    }
  }
  printf("%d file(s) removed.\n\n", removed);
}

// ─── printSimulationSummary ─────────────────────────────────────────────────

static void printSimulationSummary(const Parameters* params) {
  print_separator("═");
  printf("  SIMULATION COMPLETE\n");
  print_separator("═");

  FILE* po = std::fopen("po.dat", "r");
  if (po) {
    double sumS = 0, sumE = 0;
    double maxS = -1e9, minS = 1e9;
    double maxE = -1e9, minE = 1e9;
    int count = 0;
    char line[256];

    if (std::fgets(line, sizeof(line), po) != nullptr) {
      double T, S, varS, E, varE;
      while (std::fgets(line, sizeof(line), po) != nullptr) {
        if (std::sscanf(line, "%lf %lf %lf %lf %lf", &T, &S, &varS, &E, &varE) == 5) {
          sumS += S;  sumE += E;
          if (S > maxS) maxS = S;  if (S < minS) minS = S;
          if (E > maxE) maxE = E;  if (E < minE) minE = E;
          count++;
        }
      }
    }
    std::fclose(po);

    if (count > 0) {
      printf("\n  %-28s  mean=%8.4f   max=%8.4f   min=%8.4f\n",
             "Order parameter S:", sumS/count, maxS, minS);
      printf("  %-28s  mean=%8.4f   max=%8.4f   min=%8.4f\n",
             "Energy E:", sumE/count, maxE, minE);
      printf("  %-28s  %d\n", "Saved steps:", count);
    } else {
      printf("  [!] No valid data found in po.dat\n");
    }
  } else {
    printf("  [!] po.dat not found\n");
  }

  printf("\n");
  print_separator("─");
  printf("  %-20s %dx%dx%d\n",  "Grid:",       params->Nx, params->Ny, params->Nz);
  printf("  %-20s %s\n",        "Device:",     params->use_gpu ? "GPU (CUDA)" : "CPU (OpenMP)");
  printf("  %-20s %s\n",        "Potential:",  params->potential);
  printf("  %-20s %s\n",        "Geometry:",   params->geometry);
  printf("  %-20s %s / %s / %s\n", "Boundaries:",
         params->XBoundtype, params->YBoundtype, params->ZBoundtype);
  print_separator("═");
  printf("\n");
}

// ─── parseArgs ──────────────────────────────────────────────────────────────

static const char* parseArgs(int argc, char** argv, int* gpu_override) {
  *gpu_override = -1;
  const char* param_file = nullptr;

  for (int i = 1; i < argc; i++) {
    if (std::strcmp(argv[i], "--gpu") == 0) {
      *gpu_override = 1;
    } else if (std::strcmp(argv[i], "--cpu") == 0) {
      *gpu_override = 0;
    } else if (argv[i][0] == '-') {
      std::fprintf(stderr, "  [ERROR] Unknown option: %s\n", argv[i]);
      return nullptr;
    } else {
      if (param_file) {
        std::fprintf(stderr, "  [ERROR] Multiple parameter files specified.\n");
        return nullptr;
      }
      param_file = argv[i];
    }
  }

  return param_file;
}

// ─── main ───────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {
  print_banner();

  int gpu_override = -1;
  const char* param_file = parseArgs(argc, argv, &gpu_override);

  if (!param_file) {
    printf("  Usage: %s <parameter_file> [--gpu | --cpu]\n\n", argv[0]);
    printf("  Options:\n");
    printf("    --gpu   Force GPU execution (requires CUDA build)\n");
    printf("    --cpu   Force CPU execution (overrides 'device gpu' in param file)\n");
    printf("    Device can also be set in the parameter file: 'device gpu' or 'device cpu'\n\n");
    return 1;
  }

  cleanPreviousOutput();

  Parameters params = read_input_file(param_file);

  if (gpu_override == 1) {
    params.use_gpu = true;
    printf("  [CLI] Device override → GPU\n\n");
  } else if (gpu_override == 0) {
    params.use_gpu = false;
    printf("  [CLI] Device override → CPU\n\n");
  }

  print_parameters(params);

  simulator* sim = new simulator(&params);
  sim->Setup_simmulation(params);

  sim->print_n("ic.csv", &params);
  sim->evolve->run();

  printSimulationSummary(&params);

  delete sim;
  return 0;
}
