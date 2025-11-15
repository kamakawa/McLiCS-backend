#ifndef IC_H_
#define IC_H_
#include <gsl/gsl_eigen.h>
#include <iostream>

#include "../include/define.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"

namespace InitialConditions {
    void random(float* ni, int* pt, Parameters params);
    void homogeneous(float* ni, int* pt, Parameters params);
    void cholesteric(float* ni, int* pt, Parameters params);
    void lhelix(float* ni, int* pt, Parameters params);
    void fromFile(float* ni, int* pt, Parameters params);
    void apply(float* ni, int* pt, Parameters params);
}

// Backward compatibility - mantém exatamente as mesmas funções originais
inline void random_ic(float* ni, int* pt, Parameters params) {
    InitialConditions::random(ni, pt, params);
}

inline void homogeneous_ic(float* ni, int* pt, Parameters params) {
    InitialConditions::homogeneous(ni, pt, params);
}

inline void cholesteric_ic(float* ni, int* pt, Parameters params) {
    InitialConditions::cholesteric(ni, pt, params);
}

inline void lhelix_ic(float* ni, int* pt, Parameters params) {
    InitialConditions::lhelix(ni, pt, params);
}

inline void read_ic_file(float* ni, int* pt, Parameters params) {
    InitialConditions::fromFile(ni, pt, params);
}

inline void apply_Initial_Conditions(float* ni, int* pt, Parameters params) {
    InitialConditions::apply(ni, pt, params);
}

#endif