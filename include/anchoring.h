#ifndef ANCHORING_H_
#define ANCHORING_H_

// --- System Includes ---
#include <gsl/gsl_rng.h>
#include <iostream>
#include <string>
#include <vector>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

class Anchoring {
 public:
  // Virtual destructor ensures derived classes are properly destroyed
  virtual ~Anchoring() = default;

  char name[50];
  int id;
  Parameters* params;

  virtual char* getName() { return name; }
  
  // Base surface potential (defaults to 0)
  virtual float surface_potential(float ni[3], float s[3]) { return 0; }
  
  void check_parameter(bool std_val, std::string parameter_name);
};

class Strong_Anchoring : public Anchoring {
 public:
  char name[50] = "Strong Anchoring";
  
  float W, theta_s, phi_s;

  Strong_Anchoring(Parameters* params, int id);

  char* getName() override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
};

class RP_Anchoring : public Anchoring {
 public:
  char name[50] = "Rapine Papoular Anchoring";
  
  float W, theta_s, phi_s;

  RP_Anchoring(Parameters* params, int id);

  char* getName() override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
};

class FG_Anchoring : public Anchoring {
 public:
  char name[50] = "Founier-Galatola like Anchoring";
  
  float W, theta_s, phi_s;

  FG_Anchoring(Parameters* params, int id);

  char* getName() override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
};

class Homeotropic_Anchoring : public Anchoring {
 public:
  char name[50] = "Homeotropic Anchoring";
  
  float W, theta_s, phi_s;

  Homeotropic_Anchoring(Parameters* params, int id);

  char* getName() override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
};

class Strong_Anchoring_GHRL : public Anchoring {
 public:
  char name[50] = "Strong Anchoring GHRL";
  
  float W, theta_s, phi_s;

  Strong_Anchoring_GHRL(Parameters* params, int id);

  char* getName() override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
};

class RP_Anchoring_GHRL : public Anchoring {
 public:
  char name[50] = "Rapine Papoular Anchoring GHRL";
  
  float W, theta_s, phi_s;

  RP_Anchoring_GHRL(Parameters* params, int id);

  char* getName() override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
};

class FG_Anchoring_GHRL : public Anchoring {
 public:
  char name[50] = "Founier-Galatola like Anchoring GHRL";
  
  float W, theta_s, phi_s;

  FG_Anchoring_GHRL(Parameters* params, int id);

  char* getName() override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
};

class Homeotropic_Anchoring_GHRL : public Anchoring {
 public:
  char name[50] = "Homeotropic Anchoring GHRL";
  
  float W, theta_s, phi_s;

  Homeotropic_Anchoring_GHRL(Parameters* params, int id);

  char* getName() override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
};

#endif