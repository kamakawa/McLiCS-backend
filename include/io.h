#ifndef IO_H_
#define IO_H_

// --- System Includes ---
#include <fstream> 

// --- Project Includes ---
#include "../include/define.h"
#include "../include/parameters.h"

// Output Functions (Printing)
int print_nbc(char *fname, float *ni, float *bi, float *ci, Parameters params);

// Imprime o campo diretor principal e tipos de pontos
int print_n(char *fname, float *ni, Parameters params, int *pt);

// Imprime os parametros atuais no console/log
void print_parameters(Parameters params);

// Input & Setup Functions
// Configura parametros especificos para simulacoes GHRL
void setGHRL(Parameters &params);

// Le o arquivo de entrada de parametros
Parameters read_input_file(char *fname);

// Configura a simulacao (alocacao de memoria e estado inicial)
void Setup_simmulation(float **ni, float **bi, float **ci, Parameters &params);

// Verifica erros na leitura do arquivo de entrada
void check_error_bits(std::ifstream *f, char *parser);

#endif