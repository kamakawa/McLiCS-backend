//Seguindo o princío da POO, transformei isso em um namespace
//Por más que não envolve uma classe, segue a mesma idéia de encapsulamento e organização do código
#ifndef IO_H_
#define IO_H_
#include "../include/define.h"
#include "../include/parameters.h"

namespace IO {
    
    int print_nbc(char *fname, float *ni, float *bi, float *ci, Parameters params);
    int print_n(char *fname, float *ni, Parameters params, int *pt);

    void setGHRL(Parameters &params);
    Parameters read_input_file(char *fname);
    void print_parameters(Parameters params);
    void Setup_simmulation(float **ni, float **bi, float **ci, Parameters &params);
    void check_error_bits(std::ifstream *f, char *parser);

}

#endif