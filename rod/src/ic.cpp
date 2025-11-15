#include "../include/ic.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>

#include "../include/define.h"
#include "../include/monte_carlo.h"
#include "../include/parameters.h"

#include <memory>

using namespace std;

void apply_Initial_Condidions(std::unique_ptr<float[]>& ni, std::unique_ptr<int[]>& pt, Parameters& params) {
  if (strcasecmp(params.ic.c_str(), "random") == 0) {
    random_ic(ni, pt, params);
  } else if (strcasecmp(params.ic.c_str(), "homogeneous") == 0) {
    homogeneous_ic(ni, pt, params);
  } else if (strcasecmp(params.ic.c_str(), "ic_file") == 0) {
    read_ic_file(ni, pt, params);
  } else if (strcasecmp(params.ic.c_str(), "cholesteric") == 0) {
    cholesteric_ic(ni, pt, params);
  } else if (strcasecmp(params.ic.c_str(), "lhelix") == 0) {
    lhelix_ic(ni, pt, params);
  } else {
    fprintf(stderr, "initial condition %s not implemented \n", params.ic.c_str());
    exit(2);
  }
}

void random_ic(std::unique_ptr<float[]>& ni, std::unique_ptr<int[]>& pt, Parameters& params) {
  int Nx = params.Nx;
  int Ny = params.Ny;
  int Nz = params.Nz;
  double nx, ny;
  gsl_rng *rng;
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng, 1);
  
  printf("Random initial conditions\n\n");
  
  float *ni_ptr = ni.get();
  int *pt_ptr = pt.get();

  for (int i = 0; i < params.Nx; i++) {
    for (int j = 0; j < params.Ny; j++) {
      for (int k = 0; k < params.Nz; k++) {
        // Acesso via macros que usam o ni_ptr e pt_ptr
        if (pt_ptr[i + Nx * (j + Ny * k)]) {
          double theta = gsl_rng_uniform(rng) * 2 * M_PI;
          
          nx = cos(theta);
          ny = sin(theta);
          
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 0] = nx;
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 1] = ny;
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 2] = 0;
        } else {
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 0] = 0;
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 1] = 0;
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 2] = 1;
        }
      }
    }
  }
  
  gsl_rng_free(rng);
}

void homogeneous_ic(std::unique_ptr<float[]>& ni, std::unique_ptr<int[]>& pt, Parameters& params) {
  int Nx = params.Nx;
  int Ny = params.Ny;
  int Nz = params.Nz;

  float theta = params.theta_0;
  float phi = params.phi_0;
  printf("homogeneous initial conditions\n\n");
  printf("phi0:   %g\n", phi);
  printf("theta0: %g\n", theta);
  
  float *ni_ptr = ni.get();
  int *pt_ptr = pt.get();

  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
        if (pt_ptr[i + Nx * (j + Ny * k)]) {
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 0] = sin(theta * M_PI / 180) * cos(phi * M_PI / 180);
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 1] = sin(theta * M_PI / 180) * sin(phi * M_PI / 180);
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 2] = cos(theta * M_PI / 180);
        } else {
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 0] = 0;
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 1] = 0;
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 2] = 1;
        }
      }
    }
  }
}

void read_ic_file(std::unique_ptr<float[]>& ni, std::unique_ptr<int[]>& pt, Parameters& params) {
  int Nx = params.Nx;
  int Ny = params.Ny;
  int Nz = params.Nz;

  FILE *input_file = fopen(params.ic_file.c_str(), "r");
  char *garbage = (char *)malloc(500);
  if (input_file == 0) {
    perror(params.ic_file.c_str());
    exit(2);
  }
  int x, y, z;
  float nx, ny, nz, bx, by, bz, cx, cy, cz;
  int line = 1;
  x = fscanf(input_file, "x,y,z,nx,ny,nz");
  garbage = fgets(garbage, 500, input_file);
  printf("Reading initial conditions from %s\n", params.ic_file.c_str());
  fflush(stdout);

  float *ni_ptr = ni.get();

  for (int i = 0; i < Nx * Ny * Nz; i++) {
    if (fscanf(input_file, "%d,%d,%d,%f,%f,%f", &x, &y, &z, &nx, &ny, &nz) != 6) {
      fprintf(stderr, "%d,%d,%d,%f,%f,%f\n", x, y, z, nx, ny, nz);
      fprintf(stderr, "Problem in the line %d!\nAbborting!!!\n", ++line);
      exit(3);
    } else {
      ni_ptr[(x + Nx * (y + Ny * z)) * 3 + 0] = nx;
      ni_ptr[(x + Nx * (y + Ny * z)) * 3 + 1] = ny;
      ni_ptr[(x + Nx * (y + Ny * z)) * 3 + 2] = nz;
      garbage = fgets(garbage, 500, input_file);
      line++;
    }
  }
  fclose(input_file);
  free(garbage); 
}

void cholesteric_ic(std::unique_ptr<float[]>& ni, std::unique_ptr<int[]>& pt, Parameters& params) {
  int Nx = params.Nx;
  int Ny = params.Ny;
  int Nz = params.Nz;

  float theta = params.theta_0;
  float phi_0 = params.phi_0 * M_PI / 180;
  float p0 = params.p0_i ? params.p0_i : params.p0;
  float phi;
  printf("Cholesteric initial conditions\n\n");
  printf("p0:     %g\n", p0);
  printf("phi0:   %g\n", phi_0);
  printf("theta0: %g\n", theta);
  
  float *ni_ptr = ni.get();
  int *pt_ptr = pt.get();

  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
        if (pt_ptr[i + Nx * (j + Ny * k)]) {
          phi = phi_0 - 2 * M_PI * k / p0;

          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 0] = sin(theta * M_PI / 180) * cos(phi);
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 1] = sin(theta * M_PI / 180) * sin(phi);
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 2] = cos(theta * M_PI / 180);
        } else {
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 0] = 0;
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 1] = 0;
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 2] = 1;
        }
      }
    }
  }
}

void lhelix_ic(std::unique_ptr<float[]>& ni, std::unique_ptr<int[]>& pt, Parameters& params) {
  int Nx = params.Nx;
  int Ny = params.Ny;
  int Nz = params.Nz;

  float theta = params.theta_0;
  float phi_0 = params.phi_0 * M_PI / 180;
  float p0 = params.p0_i ? params.p0_i : params.p0;
  float phi;
  printf("Cholesteric initial conditions\n\n");
  printf("p0:     %g\n", p0);
  printf("phi0:   %g\n", phi_0);
  printf("theta0: %g\n", theta);
  
  float *ni_ptr = ni.get();
  int *pt_ptr = pt.get();

  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
        if (pt_ptr[i + Nx * (j + Ny * k)]) {
          phi = phi_0 - 2 * M_PI * i / p0;

          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 0] = cos(theta * M_PI / 180);
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 1] = sin(theta * M_PI / 180) * cos(phi);
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 2] = sin(theta * M_PI / 180) * sin(phi);
        } else {
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 0] = 0;
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 1] = 0;
          ni_ptr[(i + Nx * (j + Ny * k)) * 3 + 2] = 1;
        }
      }
    }
  }
}