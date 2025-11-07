#ifndef MONTE_CARLO_H_
#define MONTE_CARLO_H_

#include <gsl/gsl_rng.h>
#include <memory>

#include "../include/parameters.h"

void Monte_Carlo_Step_NBC(std::unique_ptr<float[]>& ni, std::unique_ptr<float[]>& bi, std::unique_ptr<float[]>& ci, float &ang_var, gsl_rng *r, Parameters& params);
void Monte_Carlo_Step_N(std::unique_ptr<float[]>& ni, float &ang_var, gsl_rng *r, Parameters& params);

#endif