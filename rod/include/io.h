#ifndef IO_H_
#define IO_H_
#include <string> 
#include <fstream> 

#include "../include/define.h"
#include "../include/parameters.h"

int print_nbc(std::string& fname, float *ni, float *bi, float *ci, Parameters& params);

int print_n(std::string& fname, float *ni, Parameters& params, int *pt);

void setGHRL(Parameters &params);

Parameters read_input_file(std::string& fname);

void print_parameters(Parameters& params);

void Setup_simmulation(float **ni, float **bi, float **ci, Parameters &params);

void check_error_bits(std::ifstream& f, std::string& parser);

#endif