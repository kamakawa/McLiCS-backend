#include "../include/io.h"

#include <math.h>
#include <string.h>

#include <fstream>
#include <iostream>

#include "../include/define.h"
#include "../include/evolve.h"
#include "../include/ic.h"
#include "../include/parameters.h"
#include "../include/potential.h"

// ─── Helpers ────────────────────────────────────────────────────────────────

static void print_separator(const char *seg = "-", int count = 60) {
  for (int i = 0; i < count; i++) fputs(seg, stdout);
  putchar('\n');
}

static void print_header(const char *title, int count = 60) {
  print_separator("-", count);
  int titleLen = (int)strlen(title);
  int pad = (count - titleLen - 2) / 2;
  if (pad < 0) pad = 0;
  printf("%*s %s %*s\n", pad, "", title, pad, "");
  print_separator("-", count);
}

// ─── print_n ────────────────────────────────────────────────────────────────

int print_n(char *fname, float *ni, Parameters params, int *pt) {
  int Nx = params.Nx;
  int Ny = params.Ny;
  int Nz = params.Nz;
  float S;
  FILE *output = fopen(fname, "w");
  if (output == 0) {
    perror(fname);
    return 1;
  }
  fprintf(output, "x,y,z,nx,ny,nz,S,pt\n");

  for (int k = 0; k < params.Nz; k++) {
    for (int j = 0; j < params.Ny; j++) {
      for (int i = 0; i < params.Nx; i++) {
        // pt == 0 means the site is outside the simulated geometry (no
        // material there) - there is no order parameter to report. NAN
        // is used instead of a numeric placeholder like 0 or 1, since
        // both of those are valid physical S values and would silently
        // masquerade as real data in downstream analysis.
        S = pti(i, j, k) ? lattice_order_parameter(ni, pt, i, j, k, params) : NAN;
        fprintf(output, "%d,%d,%d,%g,%g,%g,%g,%d\n", i, j, k,
                nix(i, j, k), niy(i, j, k), niz(i, j, k), S, pti(i, j, k));
      }
    }
  }
  printf("  Snapshot -> %s\n", fname);
  fflush(stdout);
  fclose(output);
  return 0;
}
Parameters read_input_file(char *fname) {
  Parameters input_params;
  std::ifstream input_file;
  int nn;
  if (fname == NULL) {
    printf("  [!] No parameter file specified - using default parameter values.\n\n");
    return input_params;
  } else {
    input_file.open(fname);
    if (!input_file) {
      perror(fname);
      exit(5);
    }
    printf("  Parameter file : %s\n\n", fname);
  }
  char parser[200];
  char *garbage = (char *)malloc(400);

  while (input_file >> parser) {
    if (strcasecmp(parser, "nx") == 0) {
      input_file >> input_params.Nx;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "ny") == 0) {
      input_file >> input_params.Ny;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "nz") == 0) {
      input_file >> input_params.Nz;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "A") == 0) {
      input_file >> input_params.A;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "B1") == 0) {
      input_file >> input_params.B1;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "B2") == 0) {
      input_file >> input_params.B2;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "C") == 0) {
      input_file >> input_params.C;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "Ti") == 0) {
      input_file >> input_params.Ti;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "Tf") == 0) {
      input_file >> input_params.Tf;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "dT") == 0) {
      input_file >> input_params.dT;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "MCS") == 0) {
      input_file >> input_params.MCS;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "MCT") == 0) {
      input_file >> input_params.MCT;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "ic") == 0) {
      input_file >> input_params.ic;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "ic_file") == 0) {
      input_file >> input_params.ic_file;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "xbound") == 0) {
      input_file >> input_params.XBoundtype;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "ybound") == 0) {
      input_file >> input_params.YBoundtype;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "zbound") == 0) {
      input_file >> input_params.ZBoundtype;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "potential") == 0) {
      input_file >> input_params.potential;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "evol")*strcasecmp(parser, "mode") == 0) {
      input_file >> input_params.evol;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "k11") == 0) {
      input_file >> input_params.k11;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "k22") == 0) {
      input_file >> input_params.k22;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "k33") == 0) {
      input_file >> input_params.k33;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "p0") == 0) {
      input_file >> input_params.p0;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "phi_0") == 0) {
      input_file >> input_params.phi_0;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "theta_0") == 0) {
      input_file >> input_params.theta_0;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "p0_i") == 0) {
      input_file >> input_params.p0_i;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "fn") == 0) {
      input_file >> input_params.fn;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "first_file_number") * strcasecmp(parser, "initial_output") == 0) {
      input_file >> input_params.first_file;
      input_file.getline(garbage, 400);
    } else if (parser[0] == '#') {
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "geometry") == 0) {
      input_file >> input_params.geometry;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "boundary_file") == 0) {
      input_file >> input_params.bound_file_name;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "anchoring_type") == 0) {
      input_file >> nn;
      input_file >> input_params.anchoring_type[nn];
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "W") == 0) {
      input_file >> nn;
      input_file >> input_params.W[nn];
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "phi_s") == 0) {
      input_file >> nn;
      input_file >> input_params.phi_s[nn];
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "theta_s") == 0) {
      input_file >> nn;
      input_file >> input_params.theta_s[nn];
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "nk") * strcasecmp(parser, "neighbour_kind") * strcasecmp(parser, "neighbour") == 0) {
      int nk_val;
      input_file >> nk_val;
      input_file.getline(garbage, 400);
      if (nk_val != 1) {
        fprintf(stderr,
          "\n  [ERROR] neighbour_kind = %d is not supported.\n"
          "          Only nk = 1 (nearest neighbours) is allowed.\n"
          "          Please remove or correct the 'nk' entry in your parameter file.\n\n",
          nk_val);
        exit(1);
      }
      input_params.neighbourKind = 1;
    } else if (strcasecmp(parser, "elecX")*strcasecmp(parser, "Electric_field_x")*strcasecmp(parser, "EF_x") == 0) {
      input_file >> input_params.elecX;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "elecY")*strcasecmp(parser, "Electric_field_y")*strcasecmp(parser, "EF_y")  == 0) {
      input_file >> input_params.elecY;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "elecZ")*strcasecmp(parser, "Electric_field_z")*strcasecmp(parser, "EF_z")  == 0) {
      input_file >> input_params.elecZ;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "Ei")*strcasecmp(parser, "initial_E")  == 0) {
      input_file >> input_params.elecEi;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "Ef")*strcasecmp(parser, "final_E")  == 0) {
      input_file >> input_params.elecEf;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "dE")*strcasecmp(parser, "E_variation")  == 0) {
      input_file >> input_params.elecdE;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "elecA")*strcasecmp(parser, "dielectric_anisotropy")*strcasecmp(parser, "Aniso_E")  == 0) {
      input_file >> input_params.elecA;
      input_file.getline(garbage, 400);
    } else {
      fprintf(stderr, "  [ERROR] Unknown parameter: '%s'\n", parser);
      exit(3);
    }
    check_error_bits(&input_file, parser);
  }

  return input_params;
}
// ─── print_parameters ───────────────────────────────────────────────────────

void print_parameters(Parameters params) {
  print_header("SIMULATION PARAMETERS");
  printf("  %-20s %dx%dx%d\n",  "Grid:",      params.Nx, params.Ny, params.Nz);
  printf("  %-20s %s\n",        "Geometry:",  params.geometry);
  printf("  %-20s %s\n",        "Potential:", params.potential);
  printf("  %-20s %s\n",        "Mode:",      params.evol);
  printf("  %-20s %d steps\n",  "MCS:",       params.MCS);
  printf("  %-20s %d steps\n",  "MCT:",       params.MCT);
  if (strcasecmp(params.ic, "ic_file") == 0)
    printf("  %-20s %s\n",      "IC file:",   params.ic_file);
  if (params.elecA) {
    printf("\n  Electric Field\n");
    printf("  %-20s (%g, %g, %g)\n", "Direction:",  params.elecX, params.elecY, params.elecZ);
    printf("  %-20s %g\n",           "Anisotropy:", params.elecA);
  }
  print_separator();
  printf("\n");
}
// ─── setGHRL ────────────────────────────────────────────────────────────────

void setGHRL(Parameters &params) {
  float k11 = params.k11;
  float k22 = params.k22;
  float k33 = params.k33;
  float p0  = params.p0;
  float Scale = fabs((9.0) / (k11 - 3 * k22 - k33));

  params.ghrl_sigma  = p0 ? -Scale * k22 * (2 * M_PI / p0) : 0;
  params.ghrl_lambda = (Scale / 9) * (2 * k11 - 3 * k22 + k33);
  params.ghrl_mu     = Scale * (k22 - k11);
  params.ghrl_rho    = (Scale / 9) * (k11 - k33);
  params.ghrl_nu     = (Scale / 9) * (k11 - 3 * k22 - k33);

  print_header("ELASTIC PARAMETERS");
  printf("  %-12s %g\n", "k11:",    k11);
  printf("  %-12s %g\n", "k22:",    k22);
  printf("  %-12s %g\n", "k33:",    k33);
  printf("  %-12s %g\n", "lambda:", params.ghrl_lambda);
  printf("  %-12s %g\n", "mu:",     params.ghrl_mu);
  printf("  %-12s %g\n", "rho:",    params.ghrl_rho);
  printf("  %-12s %g\n", "nu:",     params.ghrl_nu);
  printf("  %-12s %g\n", "sigma:",  params.ghrl_sigma);
  printf("  %-12s %g\n", "scale:",  Scale);
  print_separator();
  printf("\n");
}

// ─── check_error_bits ───────────────────────────────────────────────────────

void check_error_bits(std::ifstream *f, char *parser) {
  if (f->eof()) {
    // EOF after getline() is normal, not an error condition: in case there is
    // data between the last delimiter and EOF, getline() extracts it and
    // sets the eofbit.
    return;
  }
  if (f->fail()) {
    std::cerr << "  [ERROR] Invalid value for parameter: " << parser << std::endl;
    exit(1);
  }
  if (f->bad()) {
    std::cerr << "  [ERROR] Stream error while reading parameter: " << parser << std::endl;
    exit(1);
  }
}
