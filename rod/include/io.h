#ifndef IO_H_
#define IO_H_

#include <fstream>
#include <string>
#include "../include/define.h"
#include "../include/parameters.h"

namespace FileIO {
    int print_nbc(const char* fname, const float* ni, const float* bi, const float* ci, const Parameters& params);
    int print_n(const char* fname, const float* ni, const Parameters& params, const int* pt);
    Parameters read_input_file(const char* fname);
    void print_parameters(const Parameters& params);
    void check_error_bits(std::ifstream* f, const char* parser);
}

namespace SimulationSetup {
    void setGHRL(Parameters& params);
    void setup_simulation(float** ni, float** bi, float** ci, Parameters& params);
}

// Backward compatibility - mantém exatamente as mesmas funções originais
inline int print_nbc(char* fname, float* ni, float* bi, float* ci, Parameters params) {
    return FileIO::print_nbc(fname, ni, bi, ci, params);
}

inline int print_n(char* fname, float* ni, Parameters params, int* pt) {
    return FileIO::print_n(fname, ni, params, pt);
}

inline Parameters read_input_file(char* fname) {
    return FileIO::read_input_file(fname);
}

inline void print_parameters(Parameters params) {
    FileIO::print_parameters(params);
}

inline void setGHRL(Parameters& params) {
    SimulationSetup::setGHRL(params);
}

inline void Setup_simmulation(float** ni, float** bi, float** ci, Parameters& params) {
    SimulationSetup::setup_simulation(ni, bi, ci, params);
}

inline void check_error_bits(std::ifstream* f, char* parser) {
    FileIO::check_error_bits(f, parser);
}

#endif