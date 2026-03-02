#ifndef IO_H_
#define IO_H_

// --- System Includes ---
#include <fstream>

// --- Project Includes ---
#include "../include/define.h"
#include "../include/parameters.h"

// ============================================================
// MELHORIA CRÍTICA (PERFORMANCE + CUDA-READY):
// NÃO copiar Parameters por valor (contém std::map/std::string).
// Trocar para const Parameters& mantém o mesmo resultado e evita
// cópias pesadas em toda chamada.
// ============================================================

// Output Functions (Printing)
int print_nbc(char *fname, float *ni, float *bi, float *ci, const Parameters& params);

// Imprime o campo diretor principal e tipos de pontos
int print_n(char *fname, float *ni, const Parameters& params, int *pt);

// Imprime os parametros atuais no console/log
void print_parameters(const Parameters& params);

// Input & Setup Functions
// Configura parametros especificos para simulacoes GHRL
void setGHRL(Parameters &params);

// Le o arquivo de entrada de parametros
Parameters read_input_file(char *fname);

// Verifica erros na leitura do arquivo de entrada
void check_error_bits(std::ifstream *f, char *parser);

#endif