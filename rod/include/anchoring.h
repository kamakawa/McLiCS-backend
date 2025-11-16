#ifndef ANCHORING_H_
#define ANCHORING_H_

#include <gsl/gsl_rng.h>
#include <iostream>
#include <vector>
#include <string> // Necessário para std::string
#include "../include/define.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"

class Anchoring {
 public:
  // Construtor que inicializa os membros.
  Anchoring(Parameters* p, int i) : params(p), id(i) {} 
  Anchoring() = default; // Adicionado: Construtor padrão (2-C)
  virtual ~Anchoring() = default; // Adicionado: Destrutor virtual (2-B)

  std::string name; // Alterado de char[50] (1-A)
  int id;
  Parameters* params;

  // Alterado para retornar const std::string& e marcado como const (1-A, 2-A)
  virtual const std::string& getName() const { return name; }; 
  
  // Marcado como const (2-A)
  virtual float surface_potential(float ni[3], float s[3]) const { return 0; }; 
  
  void check_parameter(bool std_val, std::string parameter_name);
};

// --- Classes Derivadas (3-A e 2-A aplicados) ---

class Strong_Anchoring : public Anchoring {
 public:
  // Inicialização do nome no construtor e chamada ao construtor base (3-A)
  Strong_Anchoring(Parameters* params, int id) : Anchoring(params, id) { name = "Strong Anchoring"; } 
  float W, theta_s, phi_s;
  // Usando override e const (2-A)
  const std::string& getName() const override { return name; }; 
  float surface_potential(float ni[3], float s[3]) const override;
};

class RP_Anchoring : public Anchoring {
 public:
  RP_Anchoring(Parameters* params, int id) : Anchoring(params, id) { name = "Rapine Papoular Anchoring"; }
  float W, theta_s, phi_s;
  const std::string& getName() const override { return name; };
  float surface_potential(float ni[3], float s[3]) const override;
};

class FG_Anchoring : public Anchoring {
 public:
  FG_Anchoring(Parameters* params, int id) : Anchoring(params, id) { name = "Founier-Galatola like Anchoring"; }
  float W, theta_s, phi_s;
  const std::string& getName() const override { return name; };
  float surface_potential(float ni[3], float s[3]) const override;
};

class Homeotropic_Anchoring : public Anchoring {
 public:
  Homeotropic_Anchoring(Parameters* params, int id) : Anchoring(params, id) { name = "Homeotropic Anchoring"; }
  float W, theta_s, phi_s;
  const std::string& getName() const override { return name; };
  float surface_potential(float ni[3], float s[3]) const override;
};

class Strong_Anchoring_GHRL : public Anchoring {
 public:
  Strong_Anchoring_GHRL(Parameters* params, int id) : Anchoring(params, id) { name = "Strong Anchoring GHRL"; }
  float W, theta_s, phi_s;
  const std::string& getName() const override { return name; };
  float surface_potential(float ni[3], float s[3]) const override;
};

class RP_Anchoring_GHRL : public Anchoring {
 public:
  RP_Anchoring_GHRL(Parameters* params, int id) : Anchoring(params, id) { name = "Rapine Papoular Anchoring GHRL"; }
  float W, theta_s, phi_s;
  const std::string& getName() const override { return name; };
  float surface_potential(float ni[3], float s[3]) const override;
};

class FG_Anchoring_GHRL : public Anchoring {
 public:
  FG_Anchoring_GHRL(Parameters* params, int id) : Anchoring(params, id) { name = "Founier-Galatola like Anchoring GHRL"; }
  float W, theta_s, phi_s;
  const std::string& getName() const override { return name; };
  float surface_potential(float ni[3], float s[3]) const override;
};

class Homeotropic_Anchoring_GHRL : public Anchoring {
 public:
  Homeotropic_Anchoring_GHRL(Parameters* params, int id) : Anchoring(params, id) { name = "Homeotropic Anchoring GHRL"; }
  float W, theta_s, phi_s;
  const std::string& getName() const override { return name; };
  float surface_potential(float ni[3], float s[3]) const override;
};

#endif