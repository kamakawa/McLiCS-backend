#include "../include/io.h"

#include <math.h>
#include <string.h>
#include <fstream>
#include <iostream>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/evolve_strategy.h"
#include "../include/ic.h"
#include "../include/parameters.h"
#include "../include/potential.h"
#include "../include/parameter_order.h"  // ✅ NECESSÁRIO: OrderParam::lattice_order_parameter

// ============================================================
// Output Functions (Printing)
// ============================================================

// Imprime snapshot completo para simulações NBC (ni, bi, ci)
// Mantém a assinatura pública declarada em io.h.
int print_nbc(char *fname, float *ni, float *bi, float *ci, const Parameters& params) {
  FILE *output = fopen(fname, "w");
  if (output == 0) {
    perror(fname);
    return 1;
  }

  fprintf(output, "x,y,z,nx,ny,nz,bx,by,bz,cx,cy,cz\n");

  const int Nx = params.Nx;
  const int Ny = params.Ny;
  const int Nz = params.Nz;

  for (int k = 0; k < Nz; k++) {
    for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
        fprintf(output, "%d,%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
                i, j, k,
                nix(i, j, k), niy(i, j, k), niz(i, j, k),
                bix(i, j, k), biy(i, j, k), biz(i, j, k),
                cix(i, j, k), ciy(i, j, k), ciz(i, j, k));
      }
    }
  }

  printf("Snapshot saved in %s\n", fname);
  fflush(stdout);
  fclose(output);
  return 0;
}

// Imprime o snapshot do sistema (vetores e parametro de ordem)
int print_n(char *fname, float *ni, const Parameters& params, int *pt) {
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
        // Calcula parametro de ordem local apenas se o ponto for valido
        S = pti(i, j, k) ? OrderParam::lattice_order_parameter(ni, pt, i, j, k, params) : 1;
        
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

// Le arquivo de configuracao e preenche a struct Parameters
Parameters read_input_file(char *fname) {
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
  // Buffer para descartar o resto da linha (comentarios)
  char *garbage = (char *)malloc(400);

  while (input_file >> parser) {
    
    // --- Dimensoes ---
    if (strcasecmp(parser, "nx") == 0) {
      input_file >> input_params.Nx;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "ny") == 0) {
      input_file >> input_params.Ny;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "nz") == 0) {
      input_file >> input_params.Nz;
      input_file.getline(garbage, 400);
    } 
    
    // --- Parametros de Potencial (Lebwohl-Lasher Gen.) ---
    else if (strcasecmp(parser, "A") == 0) {
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
    } 
    
    // --- Temperatura e Monte Carlo ---
    else if (strcasecmp(parser, "Ti") == 0) {
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
    } 
    
    // --- Condicoes Iniciais e Configuracao ---
    else if (strcasecmp(parser, "ic") == 0) {
      input_file >> input_params.ic;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "ic_file") == 0) {
      input_file >> input_params.ic_file;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "potential") == 0) {
      input_file >> input_params.potential;
      input_file.getline(garbage, 400);
    } else if ((strcasecmp(parser, "evol") == 0) || (strcasecmp(parser, "mode") == 0)) {
      input_file >> input_params.evol;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "fn") == 0) {
      input_file >> input_params.fn;
      input_file.getline(garbage, 400);
    } else if ((strcasecmp(parser, "first_file_number") == 0) || (strcasecmp(parser, "initial_output") == 0)) {
      input_file >> input_params.first_file;
      input_file.getline(garbage, 400);
    } 
    
    // --- Fronteiras (Boundaries) ---
    else if (strcasecmp(parser, "xbound") == 0) {
      input_file >> input_params.XBoundtype;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "ybound") == 0) {
      input_file >> input_params.YBoundtype;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "zbound") == 0) {
      input_file >> input_params.ZBoundtype;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "geometry") == 0) {
      input_file >> input_params.geometry;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "boundary_file") == 0) {
      input_file >> input_params.bound_file_name;
      input_file.getline(garbage, 400);
    } 
    
    // --- Ancoramento ---
    else if (strcasecmp(parser, "anchoring_type") == 0) {
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
    } 
    
    // --- Constantes Elasticas (Frank) ---
    else if (strcasecmp(parser, "k11") == 0) {
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
    } else if (strcasecmp(parser, "p0_i") == 0) {
      input_file >> input_params.p0_i;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "phi_0") == 0) {
      input_file >> input_params.phi_0;
      input_file.getline(garbage, 400);
    } else if (strcasecmp(parser, "theta_0") == 0) {
      input_file >> input_params.theta_0;
      input_file.getline(garbage, 400);
    } else if ((strcasecmp(parser, "nk") == 0) || (strcasecmp(parser, "neighbour_kind") == 0) || (strcasecmp(parser, "neighbour") == 0)) {
      input_file >> input_params.neighbourKind;
      input_file.getline(garbage, 400);
    } 
    
    // --- Campo Eletrico ---
    else if ((strcasecmp(parser, "elecX") == 0) || (strcasecmp(parser, "Electric_field_x") == 0) || (strcasecmp(parser, "EF_x") == 0)) {
      input_file >> input_params.elecX;
      input_file.getline(garbage, 400);
    } else if ((strcasecmp(parser, "elecY") == 0) || (strcasecmp(parser, "Electric_field_y") == 0) || (strcasecmp(parser, "EF_y") == 0)) {
      input_file >> input_params.elecY;
      input_file.getline(garbage, 400);
    } else if ((strcasecmp(parser, "elecZ") == 0) || (strcasecmp(parser, "Electric_field_z") == 0) || (strcasecmp(parser, "EF_z") == 0)) {
      input_file >> input_params.elecZ;
      input_file.getline(garbage, 400);
    } else if ((strcasecmp(parser, "Ei") == 0) || (strcasecmp(parser, "initial_E") == 0)) {
      input_file >> input_params.elecEi;
      input_file.getline(garbage, 400);
    } else if ((strcasecmp(parser, "Ef") == 0) || (strcasecmp(parser, "final_E") == 0)) {
      input_file >> input_params.elecEf;
      input_file.getline(garbage, 400);
    } else if ((strcasecmp(parser, "dE") == 0) || (strcasecmp(parser, "E_variation") == 0)) {
      input_file >> input_params.elecdE;
      input_file.getline(garbage, 400);
    } else if ((strcasecmp(parser, "elecA") == 0) || (strcasecmp(parser, "dielectric_anisotropy") == 0) || (strcasecmp(parser, "Aniso_E") == 0)) {
      input_file >> input_params.elecA;
      input_file.getline(garbage, 400);
    } 
    
    // --- Escalas GHRL (Scales) ---
    else if ((strcasecmp(parser, "rhoScale") == 0) || (strcasecmp(parser, "rho_Scale") == 0)) {
      input_file >> input_params.rhoScale;
      input_file.getline(garbage, 400);
    } else if ((strcasecmp(parser, "lambdaScale") == 0) || (strcasecmp(parser, "lambda_Scale") == 0)) {
      input_file >> input_params.lambdaScale;
      input_file.getline(garbage, 400);
    } else if ((strcasecmp(parser, "muScale") == 0) || (strcasecmp(parser, "mu_Scale") == 0)) {
      input_file >> input_params.muScale;
      input_file.getline(garbage, 400);
    } else if ((strcasecmp(parser, "nuScale") == 0) || (strcasecmp(parser, "nu_Scale") == 0)) {
      input_file >> input_params.nuScale;
      input_file.getline(garbage, 400);
    } else if ((strcasecmp(parser, "sigmaScale") == 0) || (strcasecmp(parser, "sigma_Scale") == 0)) {
      input_file >> input_params.sigmaScale;
      input_file.getline(garbage, 400);
    } 
    
    // --- Comentarios e Erros ---
    else if (parser[0] == '#') {
      input_file.getline(garbage, 400);
    } else {
      printf("%c\n", parser[0]);
      fprintf(stderr, "invalid parameter %s!!\n", parser);
      exit(3);
    }
    check_error_bits(&input_file, parser);
  }
  
  // Libera memoria alocada
  free(garbage);

  return input_params;
}

void print_parameters(const Parameters& params) {
  printf("###Using the following parameters for this simulation.###\n");
  printf("Nx  %d\n", params.Nx);
  printf("Ny  %d\n", params.Ny);
  printf("Nz  %d\n", params.Nz);

  printf("MCS %d\n", params.MCS);
  printf("MCT %d\n", params.MCT);
  
  if (strcasecmp(params.ic, "ic_file") == 0) 
    printf("ic_file  %s\n", params.ic_file);
  
  printf("mode  %s\n", params.potential);
  
  if (params.elecA) {
    printf("###Electric Field###\n");
    printf("Electric_Field_x  %g\n", params.elecX);
    printf("Electric_Field_y  %g\n", params.elecY);
    printf("Electric_Field_z  %g\n", params.elecZ);
    printf("Dielectric_Anisotropy  %g\n", params.elecA);
  }
  printf("\n");
}

// Calcula constantes para o potencial GHRL
void setGHRL(Parameters &params) {
  float k11 = params.k11;
  float k22 = params.k22;
  float k33 = params.k33;
  float p0  = params.p0;
  float Scale = fabs((9.0) / (k11 - 3 * k22 - k33));

  params.ghrl_sigma   = p0 ? -Scale * k22 * (2 * M_PI / p0) : 0;
  params.ghrl_lambda  = (Scale / 9) * (2 * k11 - 3 * k22 + k33);
  params.ghrl_mu      = Scale * (k22 - k11);
  params.ghrl_rho     = (Scale / 9) * (k11 - k33);
  params.ghrl_nu      = (Scale / 9) * (k11 - 3 * k22 - k33);
  params.neighbourScale = (k11 + k33) / k22;

  printf("### Elastic parameters:\n");
  printf("k11 %g \nk22 %g \nk33 %g\n", k11, k22, k33);
  printf("lambda %g\n", params.ghrl_lambda);
  printf("mu     %g\n", params.ghrl_mu);
  printf("rho    %g\n", params.ghrl_rho);
  printf("nu     %g\n", params.ghrl_nu);
  printf("sigma  %g\n\n", params.ghrl_sigma);
  printf("#scaled by  %g\n\n", Scale);

  if (params.neighbourKind == 2) {
    printf("### Second Neighbours Parameters Scale");
    printf("lambda_Scale %g\n", params.lambdaScale);
    printf("mu_Scale %g\n", params.muScale);
    printf("rho_Scale %g\n", params.rhoScale);
    printf("nu_Scale %g\n", params.nuScale);
    printf("sigma_Scale %g\n", params.sigmaScale);
  }
}

void check_error_bits(std::ifstream *f, char *parser) {
  int stop = 0;
  if (f->eof()) {
    perror("stream eofbit. error state");
    // EOF after std::getline() is not the criterion to stop processing
    // data: In case there is data between the last delimiter and EOF,
    // getline() extracts it and sets the eofbit.
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