#ifndef IO_H_
#define IO_H_
#include "../include/define.h"
#include "../include/parameters.h"

int print_nbc(char *fname, float *ni, float *bi, float *ci, Parameters params);
int print_n(char *fname, float *ni, Parameters params, int *pt);

void setGHRL(Parameters &params);
Parameters read_input_file(char *fname);
void print_parameters(Parameters params);
void Setup_simmulation(float **ni, float **bi, float **ci, Parameters &params);
void check_error_bits(std::ifstream *f, char *parser);

// ─── Terminal output helpers ────────────────────────────────────────────────
// Shared by main.cpp and every setup routine (geometry/anchoring/potential/
// evolve init) so the whole CLI output follows one consistent style instead
// of each file rolling its own separator/alignment.
void print_separator(const char *seg = "-", int count = 60);
void print_header(const char *title, int count = 60, const char *seg = "-");
// Aligned "  Label: value" line — the standard row style for setup summaries.
void print_field(const char *label, const char *value);
void print_field(const char *label, double value);
void print_field(const char *label, int value);

#endif
