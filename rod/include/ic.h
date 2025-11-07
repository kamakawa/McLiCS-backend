#ifndef _IC_H_
#define _IC_H_

#include <gsl/gsl_eigen.h>
#include <memory>

#include "../include/define.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"

void random_ic(std::unique_ptr<float[]>& ni, std::unique_ptr<int[]>& pt, Parameters& params);
void homogeneous_ic(std::unique_ptr<float[]>& ni, std::unique_ptr<int[]>& pt, Parameters& params);
void cholesteric_ic(std::unique_ptr<float[]>& ni, std::unique_ptr<int[]>& pt, Parameters& params);
void lhelix_ic(std::unique_ptr<float[]>& ni, std::unique_ptr<int[]>& pt, Parameters& params);
void read_ic_file(std::unique_ptr<float[]>& ni, std::unique_ptr<int[]>& pt, Parameters& params);
void apply_Initial_Condidions(std::unique_ptr<float[]>& ni, std::unique_ptr<int[]>& pt, Parameters& params);

#endif