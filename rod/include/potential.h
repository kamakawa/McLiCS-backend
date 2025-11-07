#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include <iostream>
#include <memory>

#include "../include/define.h"
#include "../include/parameters.h"

namespace Potential {

float Bulk_Energy_Selinger_BC(float nix, float niy, float niz, float bix, float biy, float biz, std::unique_ptr<float[]>& ni, std::unique_ptr<float[]>& bi, int i, int j, int k, Parameters& params, std::unique_ptr<int[]>& pt, int[3], int nk = 1);
float nPotential(float nix, float niy, float niz, std::unique_ptr<float[]>& ni, int i, int j, int k, Parameters& params, std::unique_ptr<int[]>& pt, float bulk_pot(float nix, float niy, float niz, std::unique_ptr<float[]>& ni, int i, int j, int k, Parameters& params, std::unique_ptr<int[]>& pt, float rij[3]));
float nbcPotential(float nix, float niy, float niz, float bix, float biy, float biz, std::unique_ptr<float[]>& ni, std::unique_ptr<float[]>& bi, int i, int j, int k, Parameters& params, std::unique_ptr<int[]>& pt,
                   float bulk_pot(float, float, float, float, float, float, std::unique_ptr<float[]>&, std::unique_ptr<float[]>&, int, int, int, Parameters&, std::unique_ptr<int[]>&, int[3]));

float Bulk_Energy_Selinger_Pear(float ni[3], float nj[3], Parameters& params, float rij[3], int nk = 1);
float Bulk_Energy_Lebwohl_Lasher(float ni[3], float nj[3], Parameters& params, float rij[3], int nk = 1);
float Bulk_Energy_GHRL(float ni[3], float nj[3], Parameters& params, float rij[3], int nk = 1);
float Electric_Potential(float ni[3], Parameters& params);

}

int Periodic_Boundary(int &ii, int NN);
int Free_Boundary(int &ii, int NN);

#endif