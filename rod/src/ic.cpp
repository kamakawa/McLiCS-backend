#include "../include/ic.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>

#include "../include/define.h"
#include "../include/monte_carlo.h"
#include "../include/parameters.h"

void apply_Initial_Condidions(float *ni, int *pt, Parameters params) {
  // Alteração 1: Acesso ao tipo de IC aninhado em 'ic'
  if (strcasecmp(params.ic.ic, "random") == 0) {
    random_ic(ni, pt, params);
  } else if (strcasecmp(params.ic.ic, "homogeneous") == 0) {
    homogeneous_ic(ni, pt, params);
  } else if (strcasecmp(params.ic.ic, "ic_file") == 0) {
    read_ic_file(ni, pt, params);
  } else if (strcasecmp(params.ic.ic, "cholesteric") == 0) {
    cholesteric_ic(ni, pt, params);
  } else if (strcasecmp(params.ic.ic, "lhelix") == 0) {
    lhelix_ic(ni, pt, params);
  } else {
    // Alteração 2: Acesso ao tipo de IC aninhado em 'ic'
    fprintf(stderr, "initial condition %s not implemented \n", params.ic.ic);
    exit(2);
  }
}

void random_ic(float *ni, int *pt, Parameters params) {
  // Alteração 3, 4, 5: Acesso às dimensões da grade aninhadas em 'lattice'
  int Nx = params.lattice.Nx;
  int Ny = params.lattice.Ny;
  int Nz = params.lattice.Nz;
  double nx, ny; // Removido componente ny
  gsl_rng *rng;
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng, 1);
  
  printf("Random initial conditions\n\n");
  
  // Alteração 6, 7, 8: Acesso às dimensões da grade aninhadas em 'lattice'
  for (int i = 0; i < params.lattice.Nx; i++) {
    for (int j = 0; j < params.lattice.Ny; j++) {
      for (int k = 0; k < params.lattice.Nz; k++) {
        if (pti(i, j, k)) {
          // Gerar um ângulo aleatório
          double theta = gsl_rng_uniform(rng) * 2 * M_PI; // ângulo entre 0 e 2π
          
          // Calcula os componentes de direção com base no ângulo
          nx = cos(theta); // Componente X
          ny = sin(theta); // Componente Z
          
          // Atribuição dos valores
          nix(i, j, k) = nx;
          niy(i, j, k) = ny; // Mantemos niy como 0
          niz(i, j, k) = 0; 
        } else {
          nix(i, j, k) = 0;
          niy(i, j, k) = 0; // Mantém niy como 0
          niz(i, j, k) = 1; // Valor fixo para o eixo Z
        }
      }
    }
  }
  
  gsl_rng_free(rng);
}

void homogeneous_ic(float *ni, int *pt, Parameters params) {
  // Alteração 9, 10, 11: Acesso às dimensões da grade aninhadas em 'lattice'
  int Nx = params.lattice.Nx;
  int Ny = params.lattice.Ny;
  int Nz = params.lattice.Nz;

  // Alteração 12, 13: Acesso a theta_0 e phi_0 aninhados em 'ic'
  float theta = params.ic.theta_0;
  float phi = params.ic.phi_0;
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

void read_ic_file(float *ni, int *pt, Parameters params) {
  // Alteração 14, 15, 16: Acesso às dimensões da grade aninhadas em 'lattice'
  int Nx = params.lattice.Nx;
  int Ny = params.lattice.Ny;
  int Nz = params.lattice.Nz;
  
  // Alteração 17: Acesso a ic_file aninhado em 'ic'
  FILE *input_file = fopen(params.ic.ic_file, "r");
  char *garbage = (char *)malloc(500);
  if (input_file == 0) {
    perror(params.ic.ic_file); // Alteração 18: Acesso a ic_file aninhado em 'ic'
    exit(2);
  }
  int x, y, z;
  float nx, ny, nz, bx, by, bz, cx, cy, cz;
  int line = 1;
  x = fscanf(input_file, "x,y,z,nx,ny,nz");
  garbage = fgets(garbage, 500, input_file);
  printf("Reading initial conditions from %s\n", params.ic.ic_file); // Alteração 19: Acesso a ic_file aninhado em 'ic'
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
      garbage = fgets(garbage, 500, input_file);
      line++;
    }
  }
  fclose(input_file);
}

void cholesteric_ic(float *ni, int *pt, Parameters params) {
  // Alteração 20, 21, 22: Acesso às dimensões da grade aninhadas em 'lattice'
  int Nx = params.lattice.Nx;
  int Ny = params.lattice.Ny;
  int Nz = params.lattice.Nz;

  // Alteração 23, 24, 25: Acesso a theta_0, phi_0 e p0_i/p0 aninhados em 'ic'
  float theta = params.ic.theta_0;
  float phi_0 = params.ic.phi_0 * M_PI / 180;
  float p0 = params.ic.p0_i ? params.ic.p0_i : params.ic.p0;
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

void lhelix_ic(float *ni, int *pt, Parameters params) {
  // Alteração 26, 27, 28: Acesso às dimensões da grade aninhadas em 'lattice'
  int Nx = params.lattice.Nx;
  int Ny = params.lattice.Ny;
  int Nz = params.lattice.Nz;

  // Alteração 29, 30, 31: Acesso a theta_0, phi_0 e p0_i/p0 aninhados em 'ic'
  float theta = params.ic.theta_0;
  float phi_0 = params.ic.phi_0 * M_PI / 180;
  float p0 = params.ic.p0_i ? params.ic.p0_i : params.ic.p0;
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