#include "../include/io.h"

#include <math.h>
#include <string.h>

#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <string>

#include "../include/define.h"
#include "../include/evolve.h"
#include "../include/ic.h"
#include "../include/parameters.h"
#include "../include/potential.h"

int print_n(const char* fname, float* ni, Parameters params, int* pt) {
  int Nx = params.Nx;
  int Ny = params.Ny;
  int Nz = params.Nz;
  float S;
  FILE* output = fopen(fname, "w");
  if (!output) {
    perror(fname);
    return 1;
  }
  fprintf(output, "x,y,z,nx,ny,nz,S,pt\n");
  for (int k = 0; k < params.Nz; k++) {
    for (int j = 0; j < params.Ny; j++) {
      for (int i = 0; i < params.Nx; i++) {
        S = pti(i, j, k) ? lattice_order_parameter(ni, pt, i, j, k, params) : 1;
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

// ---------------------------------------------------------------------------
// Input file parser — dispatch table replaces the original chain of ~60
// else-if blocks. Behaviour is identical: each key maps to a lambda that
// reads the next token(s) from the stream into the corresponding field.
// ---------------------------------------------------------------------------

using Setter = std::function<void(std::ifstream&, Parameters&)>;

// Helper: read a numeric/POD scalar into a field
template<typename T>
static Setter scalar(T Parameters::*field) {
  return [field](std::ifstream& f, Parameters& p) { f >> p.*field; };
}

// Helper: read a word into a char[N] field.
// scalar<T> cannot handle char[] because pointer-to-member of array type
// is not dereferenceable via operator>>.
template<size_t N>
static Setter cstr(char (Parameters::*field)[N]) {
  return [field](std::ifstream& f, Parameters& p) {
    std::string tmp; f >> tmp;
    strncpy(p.*field, tmp.c_str(), N - 1);
    (p.*field)[N - 1] = '\0';
  };
}

// Helper: read "index value" into a std::map<int,T> field
template<typename T>
static Setter mapEntry(std::map<int,T> Parameters::*field) {
  return [field](std::ifstream& f, Parameters& p) {
    int idx; f >> idx >> (p.*field)[idx];
  };
}

// Build the dispatch table once.
// Keys are matched case-insensitively via strcasecmp.
static const std::vector<std::pair<std::vector<std::string>, Setter>> kHandlers = {
  { {"nx"},                                           scalar(&Parameters::Nx)   },
  { {"ny"},                                           scalar(&Parameters::Ny)   },
  { {"nz"},                                           scalar(&Parameters::Nz)   },
  { {"a"},                                            scalar(&Parameters::A)    },
  { {"b1"},                                           scalar(&Parameters::B1)   },
  { {"b2"},                                           scalar(&Parameters::B2)   },
  { {"c"},                                            scalar(&Parameters::C)    },
  { {"ti"},                                           scalar(&Parameters::Ti)   },
  { {"tf"},                                           scalar(&Parameters::Tf)   },
  { {"dt"},                                           scalar(&Parameters::dT)   },
  { {"mcs"},                                          scalar(&Parameters::MCS)  },
  { {"mct"},                                          scalar(&Parameters::MCT)  },
  { {"fn"},                                           scalar(&Parameters::fn)   },
  { {"first_file_number","initial_output"},           scalar(&Parameters::first_file) },
  { {"ic"},                                           cstr(&Parameters::ic)     },
  { {"ic_file"},                                      cstr(&Parameters::ic_file) },
  { {"xbound"},                                       cstr(&Parameters::XBoundtype) },
  { {"ybound"},                                       cstr(&Parameters::YBoundtype) },
  { {"zbound"},                                       cstr(&Parameters::ZBoundtype) },
  { {"potential"},                                    cstr(&Parameters::potential)  },
  { {"evol","mode"},                                  cstr(&Parameters::evol)   },
  { {"geometry"},                                     cstr(&Parameters::geometry) },
  { {"boundary_file"},                                cstr(&Parameters::bound_file_name) },
  { {"k11"},                                          scalar(&Parameters::k11)  },
  { {"k22"},                                          scalar(&Parameters::k22)  },
  { {"k33"},                                          scalar(&Parameters::k33)  },
  { {"p0"},                                           scalar(&Parameters::p0)   },
  { {"phi_0"},                                        scalar(&Parameters::phi_0) },
  { {"theta_0"},                                      scalar(&Parameters::theta_0) },
  { {"p0_i"},                                         scalar(&Parameters::p0_i) },
  { {"nk","neighbour_kind","neighbour"},              scalar(&Parameters::neighbourKind) },
  { {"rhoscale","rho_scale"},                         scalar(&Parameters::rhoScale)      },
  { {"lambdascale","lambda_scale"},                   scalar(&Parameters::lambdaScale)   },
  { {"muscale","mu_scale"},                           scalar(&Parameters::muScale)       },
  { {"nuscale","nu_scale"},                           scalar(&Parameters::nuScale)       },
  { {"sigmascale","sigma_scale"},                     scalar(&Parameters::sigmaScale)    },
  { {"elecx","electric_field_x","ef_x"},              scalar(&Parameters::elecX)  },
  { {"elecy","electric_field_y","ef_y"},              scalar(&Parameters::elecY)  },
  { {"elecz","electric_field_z","ef_z"},              scalar(&Parameters::elecZ)  },
  { {"eleca","dielectric_anisotropy","aniso_e"},      scalar(&Parameters::elecA)  },
  { {"ei","initial_e"},                               scalar(&Parameters::elecEi) },
  { {"ef","final_e"},                                 scalar(&Parameters::elecEf) },
  { {"de","e_variation"},                             scalar(&Parameters::elecdE) },
  { {"anchoring_type"},                               mapEntry(&Parameters::anchoring_type) },
  { {"w"},                                            mapEntry(&Parameters::W)   },
  { {"phi_s"},                                        mapEntry(&Parameters::phi_s) },
  { {"theta_s"},                                      mapEntry(&Parameters::theta_s) },
};

Parameters read_input_file(const char* fname) {
  Parameters p;
  if (!fname) {
    printf("No file specified — using default parameter values\n");
    return p;
  }

  std::ifstream f(fname);
  if (!f) { perror(fname); exit(5); }
  std::cout << "Using \"" << fname << "\" as input file\n\n";

  char garbage[400];
  char token[200];

  while (f >> token) {
    // Comment lines
    if (token[0] == '#') { f.getline(garbage, sizeof(garbage)); continue; }

    // Look up the token (case-insensitive)
    bool found = false;
    for (auto& [keys, setter] : kHandlers) {
      for (auto& key : keys) {
        if (strcasecmp(token, key.c_str()) == 0) {
          setter(f, p);
          f.getline(garbage, sizeof(garbage));
          found = true;
          break;
        }
      }
      if (found) break;
    }

    if (!found) {
      fprintf(stderr, "Invalid parameter: %s\n", token);
      exit(3);
    }

    check_error_bits(&f, token);
  }
  return p;
}

void print_parameters(Parameters params) {
  printf("### Simulation parameters ###\n");
  printf("Nx  %d\n", params.Nx);
  printf("Ny  %d\n", params.Ny);
  printf("Nz  %d\n", params.Nz);
  printf("MCS %d\n", params.MCS);
  printf("MCT %d\n", params.MCT);
  if (strcasecmp(params.ic, "ic_file") == 0)
    printf("ic_file  %s\n", params.ic_file);
  printf("mode  %s\n", params.potential);
  if (params.elecA) {
    printf("### Electric field ###\n");
    printf("Electric_field_x          %g\n", params.elecX);
    printf("Electric_field_y          %g\n", params.elecY);
    printf("Electric_field_z          %g\n", params.elecZ);
    printf("Dielectric_anisotropy     %g\n", params.elecA);
  }
  printf("\n");
}

void setGHRL(Parameters& params) {
  float k11 = params.k11;
  float k22 = params.k22;
  float k33 = params.k33;
  float p0  = params.p0;
  float Scale = fabs((9.0f) / (k11 - 3*k22 - k33));

  params.ghrl_sigma  = p0 ? -Scale * k22 * (2*M_PI / p0) : 0;
  params.ghrl_lambda = (Scale / 9) * (2*k11 - 3*k22 + k33);
  params.ghrl_mu     = Scale * (k22 - k11);
  params.ghrl_rho    = (Scale / 9) * (k11 - k33);
  params.ghrl_nu     = (Scale / 9) * (k11 - 3*k22 - k33);
  params.neighbourScale = (k11 + k33) / k22;

  printf("### Elastic parameters:\n");
  printf("k11 %g\nk22 %g\nk33 %g\n", k11, k22, k33);
  printf("lambda %g\nmu     %g\nrho    %g\nnu     %g\nsigma  %g\n\n",
         params.ghrl_lambda, params.ghrl_mu, params.ghrl_rho,
         params.ghrl_nu, params.ghrl_sigma);
  printf("# scaled by %g\n\n", Scale);

  if (params.neighbourKind == 2) {
    printf("### Second-neighbour parameter scales:\n");
    printf("lambda_scale %g\n", params.lambdaScale);
    printf("mu_scale     %g\n", params.muScale);
    printf("rho_scale    %g\n", params.rhoScale);
    printf("nu_scale     %g\n", params.nuScale);
    printf("sigma_scale  %g\n", params.sigmaScale);
  }
}

void check_error_bits(std::ifstream* f, const char* parser) {
  if (f->fail()) {
    std::cerr << "Invalid value for parameter: " << parser << std::endl;
    exit(1);
  }
  if (f->bad()) {
    std::cerr << "Stream error while reading parameter: " << parser << std::endl;
    exit(1);
  }
}
