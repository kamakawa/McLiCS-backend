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

#include <memory>
#include "../include/parameter_order.h"

namespace IO {

int print_n(const std::string& fname, std::unique_ptr<float[]>& ni, Parameters& params, std::unique_ptr<int[]>& pt) {
  
  float *ni_ptr = ni.get();
  int *pt_ptr = pt.get();

  int Nx = params.Nx;
  int Ny = params.Ny;
  int Nz = params.Nz;
  float S;
  
  FILE *output = fopen(fname.c_str(), "w"); 
  if (output == 0) {
    perror(fname.c_str());
    return 1;
  }
  fprintf(output, "x,y,z,nx,ny,nz,S,pt\n");

  for (int k = 0; k < params.Nz; k++) {
    for (int j = 0; j < params.Ny; j++) {
      for (int i = 0; i < params.Nx; i++) {
        S = pti(i, j, k) ? OrderParameters::lattice_order_parameter(ni_ptr, pt_ptr, i, j, k, params) : 1;
        fprintf(output, "%d,%d,%d,%g,%g,%g,%g,%d\n", i, j, k,
                nix(i, j, k), niy(i, j, k), niz(i, j, k), S, pti(i, j, k));
      }
    }
  }
  printf("Snapshot saved in %s\n", fname.c_str());
  fflush(stdout);
  fclose(output);
  return 0;
}

Parameters read_input_file(const std::string& fname) {
  Parameters input_params;
  std::ifstream input_file;
  int nn;

  if (fname.empty()) { 
    printf("File not found, using standart parameter values\n");
    return input_params;
  } else {
    input_file.open(fname); 
    if (!input_file) {
      perror(fname.c_str());
      exit(5);
    }
    std::cout << "Using \"" << fname << "\" as input file\n\n";
  }
  
  std::string parser; 
  char *garbage = (char *)malloc(400);

  while (input_file >> parser) {
    if (strcasecmp(parser.c_str(), "nx") == 0) {
      input_file >> input_params.Nx;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "ny") == 0) {
      input_file >> input_params.Ny;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "nz") == 0) {
      input_file >> input_params.Nz;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "A") == 0) {
      input_file >> input_params.A;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "B1") == 0) {
      input_file >> input_params.B1;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "B2") == 0) {
      input_file >> input_params.B2;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "C") == 0) {
      input_file >> input_params.C;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "Ti") == 0) {
      input_file >> input_params.Ti;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "Tf") == 0) {
      input_file >> input_params.Tf;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "dT") == 0) {
      input_file >> input_params.dT;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "MCS") == 0) {
      input_file >> input_params.MCS;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "MCT") == 0) {
      input_file >> input_params.MCT;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "ic") == 0) {
      input_file >> input_params.ic;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "ic_file") == 0) {
      input_file >> input_params.ic_file;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "xbound") == 0) {
      input_file >> input_params.XBoundtype;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "ybound") == 0) {
      input_file >> input_params.YBoundtype;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "zbound") == 0) {
      input_file >> input_params.ZBoundtype;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "potential") == 0) {
      input_file >> input_params.potential;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "evol")*strcasecmp(parser.c_str(), "mode") == 0) {
      input_file >> input_params.evol;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "k11") == 0) {
      input_file >> input_params.k11;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "k22") == 0) {
      input_file >> input_params.k22;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "k33") == 0) {
      input_file >> input_params.k33;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "p0") == 0) {
      input_file >> input_params.p0;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "phi_0") == 0) {
      input_file >> input_params.phi_0;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "theta_0") == 0) {
      input_file >> input_params.theta_0;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "p0_i") == 0) {
      input_file >> input_params.p0_i;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "fn") == 0) {
      input_file >> input_params.fn;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "first_file_number") * strcasecmp(parser.c_str(), "initial_output") == 0) {
      input_file >> input_params.first_file;
      input_file.getline(garbage, 400);
    } else if (parser[0] == '#') {
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "geometry") == 0) {
      input_file >> input_params.geometry;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "boundary_file") == 0) {
      input_file >> input_params.bound_file_name;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "anchoring_type") == 0) {
      input_file >> nn;
      input_file >> input_params.anchoring_type[nn];
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "W") == 0) {
      input_file >> nn;
      input_file >> input_params.W[nn];
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "phi_s") == 0) {
      input_file >> nn;
      input_file >> input_params.phi_s[nn];
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "theta_s") == 0) {
      input_file >> nn;
      input_file >> input_params.theta_s[nn];
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "nk") * strcasecmp(parser.c_str(), "neighbour_kind") * strcasecmp(parser.c_str(), "neighbour") == 0) {
      input_file >> input_params.neighbourKind;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "elecX")*strcasecmp(parser.c_str(), "Electric_field_x")*strcasecmp(parser.c_str(), "EF_x") == 0) {
      input_file >> input_params.elecX;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "elecY")*strcasecmp(parser.c_str(), "Electric_field_y")*strcasecmp(parser.c_str(), "EF_y")  == 0) {
      input_file >> input_params.elecY;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "elecZ")*strcasecmp(parser.c_str(), "Electric_field_z")*strcasecmp(parser.c_str(), "EF_z")  == 0) {
      input_file >> input_params.elecZ;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "Ei")*strcasecmp(parser.c_str(), "initial_E")  == 0) {
      input_file >> input_params.elecEi;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "Ef")*strcasecmp(parser.c_str(), "final_E")  == 0) {
      input_file >> input_params.elecEf;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "dE")*strcasecmp(parser.c_str(), "E_variation")  == 0) {
      input_file >> input_params.elecdE;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "elecA")*strcasecmp(parser.c_str(), "dielectric_anisotropy")*strcasecmp(parser.c_str(), "Aniso_E")  == 0) {
      input_file >> input_params.elecA;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "rhoScale")*strcasecmp(parser.c_str(), "rho_Scale") == 0) {
      input_file >> input_params.rhoScale;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "lambdaScale")*strcasecmp(parser.c_str(), "lambda_Scale") == 0) {
      input_file >> input_params.lambdaScale;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "muScale")*strcasecmp(parser.c_str(), "mu_Scale") == 0) {
      input_file >> input_params.muScale;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "nuScale")*strcasecmp(parser.c_str(), "nu_Scale") == 0) {
      input_file >> input_params.nuScale;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser.c_str(), "sigmaScale")*strcasecmp(parser.c_str(), "sigma_Scale") == 0) {
      input_file >> input_params.sigmaScale;
      input_file.getline(garbage, 400);
    } else {
      printf("%s\n", parser.c_str());
      fprintf(stderr, "invalid parameter %s!!\n", parser.c_str());
      exit(3);
    }
    check_error_bits(input_file, parser); 
  }
  
  free(garbage); 

  return input_params;
}

void print_parameters(Parameters params) {
  printf("###Using the following parameters for this simulation.###\n");
  printf("Nx  %d\n", params.Nx);
  printf("Ny  %d\n", params.Ny);
  printf("Nz  %d\n", params.Nz);

  printf("MCS %d\n", params.MCS);
  printf("MCT %d\n", params.MCT);
  if (strcasecmp(params.ic.c_str(), "ic_file") == 0) printf("ic_file  %s\n", params.ic_file.c_str());
  printf("mode  %s\n", params.potential.c_str());
  if (params.elecA){
  printf("###Electric Field###\n");
  printf("Electric_Field_x  %g\n", params.elecX);
  printf("Electric_Field_y  %g\n", params.elecY);
  printf("Electric_Field_z  %g\n", params.elecZ);
  printf("Dielectric_Anisotropy  %g\n", params.elecA);


  }
  printf("\n");
}
void setGHRL(Parameters &params) {
  float k11 = params.k11;
  float k22 = params.k22;
  float k33 = params.k33;
  float p0 = params.p0;
  float Scale = fabs((9.0) / (k11 - 3 * k22 - k33));

  params.ghrl_sigma = p0 ? -Scale * k22 * (2 * M_PI / p0) : 0;
  params.ghrl_lambda = (Scale / 9) * (2 * k11 - 3 * k22 + k33);
  params.ghrl_mu = Scale * (k22 - k11);
  params.ghrl_rho = (Scale / 9) * (k11 - k33);
  params.ghrl_nu = (Scale / 9) * (k11 - 3 * k22 - k33);
  params.neighbourScale = (k11 + k33) / k22;
  printf("### Elastic parameters:\n");
  printf("k11 %g \nk22 %g \nk33 %g\n", k11, k22, k33);
  printf("lambda %g\n", params.ghrl_lambda);
  printf("mu     %g\n", params.ghrl_mu);
  printf("rho    %g\n", params.ghrl_rho);
  printf("nu     %g\n", params.ghrl_nu);
  printf("sigma  %g\n\n", params.ghrl_sigma);
  printf("#scaled by  %g\n\n", Scale);
  if (params.neighbourKind==2){
    printf("### Second Neighbours Parameters Scale");
    printf("lambda_Scale %g\n",params.lambdaScale);
    printf("mu_Scale %g\n",params.muScale);
    printf("rho_Scale %g\n",params.rhoScale);
    printf("nu_Scale %g\n",params.nuScale);
    printf("sigma_Scale %g\n",params.sigmaScale);
  }
}

void check_error_bits(std::ifstream &f, const std::string& parser) {
  int stop = 0;
  if (f.eof()) { 
    perror("stream eofbit. error state");
    stop = 0;
  }
  if (f.fail()) { 
    std::cerr << "Invalid value to parameter " << parser << std::endl;
    exit(1);
  }
  if (f.bad()) { 
    std::cerr << "stream badbit. error state" << std::endl;
    exit(1);
  }
}

} 