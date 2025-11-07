#ifndef IO_H_
#define IO_H_

#include <memory>
#include <string>

#include "../include/define.h"
#include "../include/parameters.h"

namespace IO {
    
    int print_nbc(const std::string& fname, std::unique_ptr<float[]>& ni, std::unique_ptr<float[]>& bi, std::unique_ptr<float[]>& ci, Parameters& params);
    int print_n(const std::string& fname, std::unique_ptr<float[]>& ni, Parameters& params, std::unique_ptr<int[]>& pt);

    void setGHRL(Parameters &params);
    Parameters read_input_file(const std::string& fname);
    void print_parameters(Parameters params);
    void Setup_simmulation(std::unique_ptr<float[]>& ni, std::unique_ptr<float[]>& bi, std::unique_ptr<float[]>& ci, Parameters &params);
    void check_error_bits(std::ifstream &f, const std::string& parser);

}

#endif