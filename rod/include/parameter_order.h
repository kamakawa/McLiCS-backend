#ifndef PARAMETER_ORDER_H_
#define PARAMETER_ORDER_H_

#include <gsl/gsl_eigen.h>
#include <iostream>

#include "../include/define.h"
#include "../include/parameter_order.h"  
#include "../include/parameters.h"

namespace ParameterOrder {
    float Eigen_value_evaluation(float* mat, float* vec);
    void Matrice_constructor(float* ni, float* Q, int* pt, Parameters params);
    float V_Order_parameter_evaluation(float* mat_b, float* mat_c, float* vec_b, float* vec_c, float PoB);
    float VMV(float* M, float* V);
    void C_Vector_evaluation(float* vec_n, float* vec_b, float* vec_c);
    float Polarization(float* bi, Parameters params);
    float lattice_order_parameter(float* ni, int* pt, int i, int j, int k, Parameters params);
}

// Backward compatibility 
inline float Eigen_value_evaluation(float* mat, float* vec) {
    return ParameterOrder::Eigen_value_evaluation(mat, vec);
}

inline void Matrice_constructor(float* ni, float* Q, int* pt, Parameters params) {
    ParameterOrder::Matrice_constructor(ni, Q, pt, params);
}

inline float V_Order_parameter_evaluation(float* mat_b, float* mat_c, float* vec_b, float* vec_c, float PoB) {
    return ParameterOrder::V_Order_parameter_evaluation(mat_b, mat_c, vec_b, vec_c, PoB);
}

inline float VMV(float* M, float* V) {
    return ParameterOrder::VMV(M, V);
}

inline void C_Vector_evaluation(float* vec_n, float* vec_b, float* vec_c) {
    ParameterOrder::C_Vector_evaluation(vec_n, vec_b, vec_c);
}

inline float Polarization(float* bi, Parameters params) {
    return ParameterOrder::Polarization(bi, params);
}

inline float lattice_order_parameter(float* ni, int* pt, int i, int j, int k, Parameters params) {
    return ParameterOrder::lattice_order_parameter(ni, pt, i, j, k, params);
}

#endif