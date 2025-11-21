#ifndef IC_H_
#define IC_H_

// --- System Includes ---
#include <gsl/gsl_eigen.h>
#include <iostream>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"

// Inicializa o sistema com vetores aleatorios (fase isotropica/desordenada)
void random_ic(float *ni, int *pt, Parameters params);

// Inicializa o sistema perfeitamente alinhado (fase nematica ordenada)
void homogeneous_ic(float *ni, int *pt, Parameters params);

// Inicializa o sistema com estrutura helicoidal (fase colesterica)
void cholesteric_ic(float *ni, int *pt, Parameters params);

// Le a condicao inicial de um arquivo binario ou texto
void read_ic_file(float *ni, int *pt, Parameters params);

// Funcao principal que seleciona e aplica a condicao inicial baseada nos parametros
void apply_Initial_Condidions(float *ni, int *pt, Parameters params);

#endif