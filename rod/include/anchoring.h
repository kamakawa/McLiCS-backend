#ifndef ANCHORING_H_
#define ANCHORING_H_

#include <gsl/gsl_rng.h>
#include <iostream>
#include <string>
#include <vector>
#include "../include/define.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

class Anchoring {
 public:
  virtual ~Anchoring() = default;
  virtual std::string getName() const = 0;
  virtual float surface_potential(float ni[3], float s[3]) = 0;
  void check_parameter(bool std_val, std::string parameter_name);
  
 protected:
  int getId() const { return id; }
  Parameters* getParams() const { return params; }
  
 private:
  int id;
  Parameters* params;
};

class Strong_Anchoring : public Anchoring {
 public:
  Strong_Anchoring(Parameters* params, int id);
  std::string getName() const override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
  
  // Getters para acesso controlado
  float getW() const { return W; }
  float getThetaS() const { return theta_s; }
  float getPhiS() const { return phi_s; }
  
 private:
  std::string name = "Strong Anchoring";
  float W, theta_s, phi_s;
};

class RP_Anchoring : public Anchoring {
 public:
  RP_Anchoring(Parameters* params, int id);
  std::string getName() const override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
  
  float getW() const { return W; }
  float getThetaS() const { return theta_s; }
  float getPhiS() const { return phi_s; }
  
 private:
  std::string name = "Rapine Papoular Anchoring";
  float W, theta_s, phi_s;
};

class FG_Anchoring : public Anchoring {
 public:
  FG_Anchoring(Parameters* params, int id);
  std::string getName() const override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
  
  float getW() const { return W; }
  float getThetaS() const { return theta_s; }
  float getPhiS() const { return phi_s; }
  
 private:
  std::string name = "Founier-Galatola like Anchoring";
  float W, theta_s, phi_s;
};

class Homeotropic_Anchoring : public Anchoring {
 public:
  Homeotropic_Anchoring(Parameters* params, int id);
  std::string getName() const override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
  
  float getW() const { return W; }
  float getThetaS() const { return theta_s; }
  float getPhiS() const { return phi_s; }
  
 private:
  std::string name = "Homeotropic Anchoring";
  float W, theta_s, phi_s;
};

class Strong_Anchoring_GHRL : public Anchoring {
 public:
  Strong_Anchoring_GHRL(Parameters* params, int id);
  std::string getName() const override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
  
  float getW() const { return W; }
  float getThetaS() const { return theta_s; }
  float getPhiS() const { return phi_s; }
  
 private:
  std::string name = "Strong Anchoring GHRL";
  float W, theta_s, phi_s;
};

class RP_Anchoring_GHRL : public Anchoring {
 public:
  RP_Anchoring_GHRL(Parameters* params, int id);
  std::string getName() const override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
  
  float getW() const { return W; }
  float getThetaS() const { return theta_s; }
  float getPhiS() const { return phi_s; }
  
 private:
  std::string name = "Rapine Papoular Anchoring GHRL";
  float W, theta_s, phi_s;
};

class FG_Anchoring_GHRL : public Anchoring {
 public:
  FG_Anchoring_GHRL(Parameters* params, int id);
  std::string getName() const override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
  
  float getW() const { return W; }
  float getThetaS() const { return theta_s; }
  float getPhiS() const { return phi_s; }
  
 private:
  std::string name = "Founier-Galatola like Anchoring GHRL";
  float W, theta_s, phi_s;
};

class Homeotropic_Anchoring_GHRL : public Anchoring {
 public:
  Homeotropic_Anchoring_GHRL(Parameters* params, int id);
  std::string getName() const override { return name; }
  float surface_potential(float ni[3], float s[3]) override;
  
  float getW() const { return W; }
  float getThetaS() const { return theta_s; }
  float getPhiS() const { return phi_s; }
  
 private:
  std::string name = "Homeotropic Anchoring GHRL";
  float W, theta_s, phi_s;
};

#endif