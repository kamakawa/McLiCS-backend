#ifndef POTENTIAL_H_
#define POTENTIAL_H_
#include <iostream>

#include "../include/define.h"
#include "../include/parameters.h"

namespace Potential {

    float Bulk_Energy_Selinger_BC(float nix, float niy, float niz, float bix, float biy, float biz, 
                                 float *ni, float *bi, int i, int j, int k, const Parameters& params, 
                                 int *pt, int[3], int nk = 1);
    
    float Bulk_Energy_Selinger_Pear(const float ni[3], const float nj[3], const Parameters* params, 
                                   const float rij[3], int nk = 1);
    
    float Bulk_Energy_Lebwohl_Lasher(const float ni[3], const float nj[3], const Parameters* params, 
                                    const float rij[3], int nk = 1);
    
    float Bulk_Energy_GHRL(const float ni[3], const float nj[3], const Parameters* params, 
                          const float rij[3], int nk = 1);

    float Electric_Potential(const float ni[3], const Parameters* params);
    
    float nPotential(float nix, float niy, float niz, float *ni, int i, int j, int k, 
                    const Parameters& params, int *pt, 
                    float bulk_pot(float nix, float niy, float niz, float *ni, int i, int j, int k, 
                                   const Parameters& params, int *pt, float rij[3]));
    
    float nbcPotential(float nix, float niy, float niz, float bix, float biy, float biz, 
                      float *ni, float *bi, int i, int j, int k, const Parameters& params, int *pt,
                      float bulk_pot(float, float, float, float, float, float, float *, float *, 
                                     int, int, int, const Parameters&, int *, int[3]));

    int Periodic_Boundary(int &ii, int NN);
    int Free_Boundary(int &ii, int NN);
}

#endif