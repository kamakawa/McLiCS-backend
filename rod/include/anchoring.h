#ifndef ANCHORING_H_
#define ANCHORING_H_

#include <gsl/gsl_rng.h>
#include <iostream>
#include <string>
#include "../include/define.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

class Anchoring {
 public:
  virtual ~Anchoring() = default;
  int id;
  Parameters* params;
  virtual std::string getName() { return "Undefined Anchoring"; };
  virtual float surface_potential(float ni[3], float s[3]) { return 0; };
  void check_parameter(bool std_val, std::string parameter_name);
};

class Strong_Anchoring : public Anchoring {
 public:
  Strong_Anchoring(Parameters* params, int id);
  float W, theta_s, phi_s;
  std::string getName() override { return "Strong Anchoring"; };
  float surface_potential(float ni[3], float s[3]) override;
};

class RP_Anchoring : public Anchoring {
 public:
  RP_Anchoring(Parameters* params, int id);
  float W, theta_s, phi_s;
  std::string getName() override { return "Rapine Papoular Anchoring"; };
  float surface_potential(float ni[3], float s[3]) override;
};

class FG_Anchoring : public Anchoring {
 public:
  FG_Anchoring(Parameters* params, int id);
  float W, theta_s, phi_s;
  std::string getName() override { return "Founier-Galatola like Anchoring"; };
  float surface_potential(float ni[3], float s[3]) override;
};

class Homeotropic_Anchoring : public Anchoring {
 public:
  Homeotropic_Anchoring(Parameters* params, int id);
  float W, theta_s, phi_s;
  std::string getName() override { return "Homeotropic Anchoring"; };
  float surface_potential(float ni[3], float s[3]) override;
};

class Strong_Anchoring_GHRL : public Anchoring {
 public:
  Strong_Anchoring_GHRL(Parameters* params, int id);
  float W, theta_s, phi_s;
  std::string getName() override { return "Strong Anchoring GHRL"; };
  float surface_potential(float ni[3], float s[3]) override;
};

class RP_Anchoring_GHRL : public Anchoring {
 public:
  RP_Anchoring_GHRL(Parameters* params, int id);
  float W, theta_s, phi_s;
  std::string getName() override { return "Rapine Papoular Anchoring GHRL"; };
  float surface_potential(float ni[3], float s[3]) override;
};

class FG_Anchoring_GHRL : public Anchoring {
 public:
  FG_Anchoring_GHRL(Parameters* params, int id);
  float W, theta_s, phi_s;
  std::string getName() override { return "Founier-Galatola like Anchoring GHRL"; };
  float surface_potential(float ni[3], float s[3]) override;
};

class Homeotropic_Anchoring_GHRL : public Anchoring {
 public:
  Homeotropic_Anchoring_GHRL(Parameters* params, int id);
  float W, theta_s, phi_s;
  std::string getName() override { return "Homeotropic Anchoring GHRL"; };
  float surface_potential(float ni[3], float s[3]) override;
};

#endif