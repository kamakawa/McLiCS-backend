#ifndef MONTE_CARLO_H_
#define MONTE_CARLO_H_

// --- System Includes ---
#include <gsl/gsl_rng.h>

// --- Project Includes ---
#include "../include/parameters.h"

// ============================================================
// MELHORIA CRÍTICA:
// NÃO copiar Parameters por valor (contém std::map/std::string).
// Isso é caro e atrapalha a migração para CUDA.
// Passar por referência constante mantém o mesmo resultado.
// ============================================================

// Passo de Monte Carlo para multiplos campos vetoriais (ni, bi, ci)
void Monte_Carlo_Step_NBC(float *ni, float *bi, float *ci, float &ang_var, gsl_rng *r, const Parameters& params);

// Passo de Monte Carlo padrao apenas para o campo diretor nematico (ni)
void Monte_Carlo_Step_N(float *ni, float &ang_var, gsl_rng *r, const Parameters& params);

#endif