#ifndef PARAMETER_ORDER_H_
#define PARAMETER_ORDER_H_

// --- System Includes ---
#include <gsl/gsl_eigen.h>
#include <iostream>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/parameters.h"

// Namespace para isolar funcoes de calculo de parametros de ordem
namespace OrderParam {

    // Constroi a matriz Q tensor para o calculo do parametro de ordem
    void Matrice_constructor(float *ni, float *Q, int *pt, Parameters params);

    // Calcula autovalores de uma matriz
    float Eigen_value_evaluation(float *mat, float *vec);

    // Multiplicacao Vetor-Matriz-Vetor (VMV)
    float VMV(float *M, float *V);

    // Avaliacao do vetor C (auxiliar para fases biaxiais)
    void C_Vector_evaluation(float *vec_n, float *vec_b, float *vec_c);

    // Avalia o parametro de ordem V
    float V_Order_parameter_evaluation(float *mat_b, float *mat_c, float *vec_b, float *vec_c, float PoB);

    // Calcula a polarizacao do sistema
    float Polarization(float *bi, Parameters params);

    // Calcula o parametro de ordem local em um ponto especifico (i,j,k)
    float lattice_order_parameter(float *ni, int *pt, int i, int j, int k, Parameters params);

} 

#endif