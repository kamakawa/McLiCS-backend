#ifndef POTENTIAL_H_
#define POTENTIAL_H_
#include <iostream>
#include <functional> // Necessário para std::function

#include "../include/define.h"
#include "../include/parameters.h"
#include "../include/potential.h"

namespace Potential {

    using BulkPotFunc = std::function<float(float nix, float niy, float niz, float *ni, int i, int j, int k, Parameters& params, int *pt, float rij[3])>;
    
    float Bulk_Energy_Selinger_BC(float nix, float niy, float niz, float bix, float biy, float biz, float *ni, float *bi, int i, int j, int k, Parameters& params, int *pt, int arr[3], int nk = 1);
    
    float nPotential(float nix, float niy, float niz, float *ni, int i, int j, int k, Parameters& params, int *pt, BulkPotFunc bulk_pot);
    
    int Periodic_Boundary(int &ii, int NN);
    int Free_Boundary(int &ii, int NN);
    
    using NBCBulkPotFunc = std::function<float(float, float, float, float, float, float, float *, float *, int, int, int, Parameters&, int *, int[3])>;

    float nbcPotential(float nix, float niy, float niz, float bix, float biy, float biz, float *ni, float *bi, int i, int j, int k, Parameters& params, int *pt,
                       NBCBulkPotFunc bulk_pot);

    float Bulk_Energy_Selinger_Pear(float ni[3], float nj[3], Parameters *params, float rij[3], int nk = 1);
    float Bulk_Energy_Lebwohl_Lasher(float ni[3], float nj[3], Parameters *params, float rij[3], int nk = 1);
    float Bulk_Energy_GHRL(float ni[3], float nj[3], Parameters *params, float rij[3], int nk = 1);
    float Electric_Potential(float ni[3], Parameters *params);

} 

#endif