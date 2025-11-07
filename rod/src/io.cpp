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

int IO::print_n(char *fname, float *ni, Parameters params, int *pt) {
  int Nx = params.lattice.Nx;
  int Ny = params.lattice.Ny;
  int Nz = params.lattice.Nz;
  float S;
  FILE *output = fopen(fname, "w");
  if (output == 0) {
    perror(fname);
    return 1;
  }
  fprintf(output, "x,y,z,nx,ny,nz,S,pt\n");

  for (int k = 0; k < params.lattice.Nz; k++) {
    for (int j = 0; j < params.lattice.Ny; j++) {
      for (int i = 0; i < params.lattice.Nx; i++) {
        S = pti(i, j, k) ? OrderParameters::lattice_order_parameter(ni, pt, i, j, k, params) : 1;
        fprintf(output, "%d,%d,%d,%g,%g,%g,%g,%d\n", i, j, k,
                nix(i, j, k), niy(i, j, k), niz(i, j, k), S, pti(i, j, k));
      }
    }
  }
  printf("Snapshot saved in %s\n", fname);
  fflush(stdout);
  fclose(output);
  return 0;
}

Parameters IO::read_input_file(char *fname) {
  Parameters input_params;
  std::ifstream input_file;
  int nn;
  if (fname == NULL) {
    printf("File not found, using standart parameter values\n");
    return input_params;
  } else {
    input_file.open(fname);
    if (!input_file) {
      perror(fname);
      exit(5);
    }
    std::cout << "Using \"" << fname << "\" as input file\n\n";
  }
  char parser[200];
  char *garbage = (char *)malloc(400);

  while (input_file >> parser) {
    if (strcasecmp(parser, "nx") == 0) {
      input_file >> input_params.lattice.Nx;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "ny") == 0) {
      input_file >> input_params.lattice.Ny;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "nz") == 0) {
      input_file >> input_params.lattice.Nz;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "A") == 0) {
      input_file >> input_params.potential.A;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "B1") == 0) {
      input_file >> input_params.potential.B1;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "B2") == 0) {
      input_file >> input_params.potential.B2;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "C") == 0) {
      input_file >> input_params.potential.C;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "Ti") == 0) {
      input_file >> input_params.potential.Ti;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "Tf") == 0) {
      input_file >> input_params.potential.Tf;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "dT") == 0) {
      input_file >> input_params.potential.dT;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "MCS") == 0) {
      input_file >> input_params.mc.MCS;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "MCT") == 0) {
      input_file >> input_params.mc.MCT;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "ic") == 0) {
      input_file >> input_params.ic.ic;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "ic_file") == 0) {
      input_file >> input_params.ic.ic_file;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "xbound") == 0) {
      input_file >> input_params.lattice.XBoundtype;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "ybound") == 0) {
      input_file >> input_params.lattice.YBoundtype;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "zbound") == 0) {
      input_file >> input_params.lattice.ZBoundtype;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "potential") == 0) {
      input_file >> input_params.potential.potential;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "evol") == 0 || strcasecmp(parser, "mode") == 0) {
      input_file >> input_params.mc.evol;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "k11") == 0) {
      input_file >> input_params.ic.k11;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "k22") == 0) {
      input_file >> input_params.ic.k22;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "k33") == 0) {
      input_file >> input_params.ic.k33;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "p0") == 0) {
      input_file >> input_params.ic.p0;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "phi_0") == 0) {
      input_file >> input_params.ic.phi_0;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "theta_0") == 0) {
      input_file >> input_params.ic.theta_0;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "p0_i") == 0) {
      input_file >> input_params.ic.p0_i;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "fn") == 0) {
      input_file >> input_params.mc.fn;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "first_file_number") == 0 || strcasecmp(parser, "initial_output") == 0) {
      input_file >> input_params.mc.first_file;
      input_file.getline(garbage, 400);
    } else if (parser[0] == '#') {
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "geometry") == 0) {
      input_file >> input_params.lattice.geometry;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "boundary_file") == 0) {
      input_file >> input_params.surface.bound_file_name;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "anchoring_type") == 0) {
      input_file >> nn;
      input_file >> input_params.surface.anchoring_type[nn];
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "W") == 0) {
      input_file >> nn;
      input_file >> input_params.surface.W[nn];
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "phi_s") == 0) {
      input_file >> nn;
      input_file >> input_params.surface.phi_s[nn];
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "theta_s") == 0) {
      input_file >> nn;
      input_file >> input_params.surface.theta_s[nn];
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "nk") == 0 || strcasecmp(parser, "neighbour_kind") == 0 || strcasecmp(parser, "neighbour") == 0) {
      input_file >> input_params.neighbourhood.neighbourKind;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "elecX") == 0 || strcasecmp(parser, "Electric_field_x") == 0 || strcasecmp(parser, "EF_x") == 0) {
      input_file >> input_params.electric.elecX;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "elecY") == 0 || strcasecmp(parser, "Electric_field_y") == 0 || strcasecmp(parser, "EF_y") == 0) {
      input_file >> input_params.electric.elecY;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "elecZ") == 0 || strcasecmp(parser, "Electric_field_z") == 0 || strcasecmp(parser, "EF_z") == 0) {
      input_file >> input_params.electric.elecZ;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "Ei") == 0 || strcasecmp(parser, "initial_E") == 0) {
      input_file >> input_params.electric.elecEi;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "Ef") == 0 || strcasecmp(parser, "final_E") == 0) {
      input_file >> input_params.electric.elecEf;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "dE") == 0 || strcasecmp(parser, "E_variation") == 0) {
      input_file >> input_params.electric.elecdE;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "elecA") == 0 || strcasecmp(parser, "dielectric_anisotropy") == 0 || strcasecmp(parser, "Aniso_E") == 0) {
      input_file >> input_params.electric.elecA;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "rhoScale") == 0 || strcasecmp(parser, "rho_Scale") == 0) {
      input_file >> input_params.neighbourhood.rhoScale;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "lambdaScale") == 0 || strcasecmp(parser, "lambda_Scale") == 0) {
      input_file >> input_params.neighbourhood.lambdaScale;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "muScale") == 0 || strcasecmp(parser, "mu_Scale") == 0) {
      input_file >> input_params.neighbourhood.muScale;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "nuScale") == 0 || strcasecmp(parser, "nu_Scale") == 0) {
      input_file >> input_params.neighbourhood.nuScale;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "sigmaScale") == 0 || strcasecmp(parser, "sigma_Scale") == 0) {
      input_file >> input_params.neighbourhood.sigmaScale;
      input_file.getline(garbage, 400);
    } else {
      printf("%c\n", parser[0]);
      fprintf(stderr, "invalid parameter %s!!\n", parser);
      exit(3);
    }
    IO::check_error_bits(&input_file, parser);
  }

  return input_params;
}

void IO::print_parameters(Parameters params) {
  printf("###Using the following parameters for this simulation.###\n");
  printf("Nx  %d\n", params.lattice.Nx);
  printf("Ny  %d\n", params.lattice.Ny);
  printf("Nz  %d\n", params.lattice.Nz);

  printf("MCS %d\n", params.mc.MCS);
  printf("MCT %d\n", params.mc.MCT);
  
  if (strcasecmp(params.ic.ic_file, "ic_file") == 0) printf("ic_file  %s\n", params.ic.ic_file);
  printf("mode  %s\n", params.potential.potential);
  
  if (params.electric.elecA){
  printf("###Electric Field###\n");
  printf("Electric_Field_x  %g\n", params.electric.elecX);
  printf("Electric_Field_y  %g\n", params.electric.elecY);
  printf("Electric_Field_z  %g\n", params.electric.elecZ);
  printf("Dielectric_Anisotropy  %g\n", params.electric.elecA);
  }
  printf("\n");
}

void IO::setGHRL(Parameters &params) {
  float k11 = params.ic.k11;
  float k22 = params.ic.k22;
  float k33 = params.ic.k33;
  float p0 = params.ic.p0;
  float Scale = fabs((9.0) / (k11 - 3 * k22 - k33));

  params.potential.ghrl_sigma = p0 ? -Scale * k22 * (2 * M_PI / p0) : 0;
  params.potential.ghrl_lambda = (Scale / 9) * (2 * k11 - 3 * k22 + k33);
  params.potential.ghrl_mu = Scale * (k22 - k11);
  params.potential.ghrl_rho = (Scale / 9) * (k11 - k33);
  params.potential.ghrl_nu = (Scale / 9) * (k11 - 3 * k22 - k33);
  
  params.neighbourhood.neighbourScale = (k11 + k33) / k22;
  printf("### Elastic parameters:\n");
  printf("k11 %g \nk22 %g \nk33 %g\n", k11, k22, k33);
  printf("lambda %g\n", params.potential.ghrl_lambda);
  printf("mu     %g\n", params.potential.ghrl_mu);
  printf("rho    %g\n", params.potential.ghrl_rho);
  printf("nu     %g\n", params.potential.ghrl_nu);
  printf("sigma  %g\n\n", params.potential.ghrl_sigma);
  printf("#scaled by  %g\n\n", Scale);
  
  if (params.neighbourhood.neighbourKind==2){
    printf("### Second Neighbours Parameters Scale\n");
    printf("lambda_Scale %g\n",params.neighbourhood.lambdaScale);
    printf("mu_Scale %g\n",params.neighbourhood.muScale);
    printf("rho_Scale %g\n",params.neighbourhood.rhoScale);
    printf("nu_Scale %g\n",params.neighbourhood.nuScale);
    printf("sigma_Scale %g\n",params.neighbourhood.sigmaScale);
  }
}

void IO::check_error_bits(std::ifstream *f, char *parser) {
  int stop = 0;
  if (f->eof()) {
    perror("stream eofbit. error state");
    stop = 0;
  }
  if (f->fail()) {
    std::cerr << "Invalid value to parameter " << parser << std::endl;
    exit(1);
  }
  if (f->bad()) {
    std::cerr << "stream badbit. error state" << std::endl;
    exit(1);
  }
}