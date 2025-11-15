#include "../include/ic.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>

#include "../include/define.h"
#include "../include/monte_carlo.h"
#include "../include/parameters.h"

// Usando o namespace para compatibilidade com header refatorado
namespace InitialConditions {

void apply(float *ni, int *pt, Parameters params) {
  if (strcasecmp(params.ic, "random") == 0) {
    random(ni, pt, params);
  } else if (strcasecmp(params.ic, "homogeneous") == 0) {
    homogeneous(ni, pt, params);
  } else if (strcasecmp(params.ic, "ic_file") == 0) {
    fromFile(ni, pt, params);
  } else if (strcasecmp(params.ic, "cholesteric") == 0) {
    cholesteric(ni, pt, params);
  } else if (strcasecmp(params.ic, "lhelix") == 0) {
    lhelix(ni, pt, params);
  } else {
    fprintf(stderr, "initial condition %s not implemented \n", params.ic);
    exit(2);
  }
}

void random(float *ni, int *pt, Parameters params) {
  int Nx = params.Nx;
  int Ny = params.Ny;
  int Nz = params.Nz;
  double nx, ny;
  gsl_rng *rng;
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng, 1);
  
  printf("Random initial conditions\n\n");
  
  for (int i = 0; i < params.Nx; i++) {
    for (int j = 0; j < params.Ny; j++) {
      for (int k = 0; k < params.Nz; k++) {
        if (pti(i, j, k)) {
          double theta = gsl_rng_uniform(rng) * 2 * M_PI;
          
          nx = cos(theta);
          ny = sin(theta);
          
          nix(i, j, k) = nx;
          niy(i, j, k) = ny;
          niz(i, j, k) = 0; 
        } else {
          nix(i, j, k) = 0;
          niy(i, j, k) = 0;
          niz(i, j, k) = 1;
        }
      }
    }
  }
  
  gsl_rng_free(rng);
}

void homogeneous(float *ni, int *pt, Parameters params) {
  int Nx = params.Nx;
  int Ny = params.Ny;
  int Nz = params.Nz;

  float theta = params.theta_0;
  float phi = params.phi_0;
  printf("homogeneous initial conditions\n\n");
  printf("phi0:   %g\n", phi);
  printf("theta0: %g\n", theta);
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
        if (pti(i, j, k)) {
          nix(i, j, k) = sin(theta * M_PI / 180) * cos(phi * M_PI / 180);
          niy(i, j, k) = sin(theta * M_PI / 180) * sin(phi * M_PI / 180);
          niz(i, j, k) = cos(theta * M_PI / 180);
        } else {
          nix(i, j, k) = 0;
          niy(i, j, k) = 0;
          niz(i, j, k) = 1;
        }
      }
    }
  }
}

void fromFile(float *ni, int *pt, Parameters params) {
  int Nx = params.Nx;
  int Ny = params.Ny;
  int Nz = params.Nz;
  FILE *input_file = fopen(params.ic_file, "r");
  char garbage[500];  // CORREÇÃO: array na stack em vez de malloc
  if (input_file == 0) {
    perror(params.ic_file);
    exit(2);
  }
  int x, y, z;
  float nx, ny, nz;
  int line = 1;
  
  // Lê o cabeçalho
  int header_read = fscanf(input_file, "x,y,z,nx,ny,nz");
  fgets(garbage, 500, input_file);  // Lê o resto da linha
  
  printf("Reading initial conditions from %s\n", params.ic_file);
  fflush(stdout);
  
  for (int i = 0; i < Nx * Ny * Nz; i++) {
    if (fscanf(input_file, "%d,%d,%d,%f,%f,%f", &x, &y, &z, &nx, &ny, &nz) != 6) {
      fprintf(stderr, "%d,%d,%d,%f,%f,%f\n", x, y, z, nx, ny, nz);
      fprintf(stderr, "Problem in the line %d!\nAbborting!!!\n", ++line);
      exit(3);
    } else {
      nix(x, y, z) = nx;
      niy(x, y, z) = ny;
      niz(x, y, z) = nz;
      fgets(garbage, 500, input_file);  // Lê o resto da linha
      line++;
    }
  }
  fclose(input_file);
  // CORREÇÃO: Não precisa de free - liberação automática
}

void cholesteric(float *ni, int *pt, Parameters params) {
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
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
        if (pti(i, j, k)) {
          phi = phi_0 - 2 * M_PI * k / p0;

          nix(i, j, k) = sin(theta * M_PI / 180) * cos(phi);
          niy(i, j, k) = sin(theta * M_PI / 180) * sin(phi);
          niz(i, j, k) = cos(theta * M_PI / 180);
        } else {
          nix(i, j, k) = 0;
          niy(i, j, k) = 0;
          niz(i, j, k) = 1;
        }
      }
    }
  }
}

void lhelix(float *ni, int *pt, Parameters params) {
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
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
        if (pti(i, j, k)) {
          phi = phi_0 - 2 * M_PI * i / p0;

          nix(i, j, k) = cos(theta * M_PI / 180);
          niy(i, j, k) = sin(theta * M_PI / 180) * cos(phi);
          niz(i, j, k) = sin(theta * M_PI / 180) * sin(phi);
        } else {
          nix(i, j, k) = 0;
          niy(i, j, k) = 0;
          niz(i, j, k) = 1;
        }
      }
    }
  }
}

} // namespace InitialConditions

// Backward compatibility - mantém as funções originais
void apply_Initial_Condidions(float *ni, int *pt, Parameters params) {
  InitialConditions::apply(ni, pt, params);
}

void random_ic(float *ni, int *pt, Parameters params) {
  InitialConditions::random(ni, pt, params);
}

void homogeneous_ic(float *ni, int *pt, Parameters params) {
  InitialConditions::homogeneous(ni, pt, params);
}

void read_ic_file(float *ni, int *pt, Parameters params) {
  InitialConditions::fromFile(ni, pt, params);
}

void cholesteric_ic(float *ni, int *pt, Parameters params) {
  InitialConditions::cholesteric(ni, pt, params);
}

void lhelix_ic(float *ni, int *pt, Parameters params) {
  InitialConditions::lhelix(ni, pt, params);
}