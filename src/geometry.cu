#include <iostream>
#include <string>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <vector>
#include "../include/anchoring.h"
#include "../include/geometry.h"
#include "../include/parameters.h"
#include "../include/potential.cuh"

#include <map>
__device__ extern float (*bulk_potentialGPU)(float ni[3], float nj[3], Parameters *params, float rij[3], int nk);
__device__ extern p_Surface_Potential *d_Surface_Potential_GPU; 
__device__ extern float *d_W; 
__device__ inline float newman_neighbours_GPU(const nni fullni[], Parameters *params){
  float rij[3];
  double E=0;
  float ni[3]={fullni[0].x,fullni[0].y,fullni[0].z};
  
  if (fullni[1].pt){
    float nj[3]={fullni[1].x,fullni[1].y,fullni[1].z};
    rij[0]=1; rij[1]=0; rij[2]=0;
    E+= bulk_potentialGPU(ni,nj,params,rij,1);
  }
  if (fullni[2].pt){
    float nj[3]={fullni[2].x,fullni[2].y,fullni[2].z};
    rij[0]=-1; rij[1]=0; rij[2]=0;
    E+= bulk_potentialGPU(ni,nj,params,rij,1);
  }
  if (fullni[3].pt){
    float nj[3]={fullni[3].x,fullni[3].y,fullni[3].z};
    rij[0]=0; rij[1]=1; rij[2]=0;
    E+= bulk_potentialGPU(ni,nj,params,rij,1);
  }
  if (fullni[4].pt){
    float nj[3]={fullni[4].x,fullni[4].y,fullni[4].z};
    rij[0]=0; rij[1]=-1; rij[2]=0;
    E+= bulk_potentialGPU(ni,nj,params,rij,1);
  }
  if (fullni[5].pt){
    float nj[3]={fullni[5].x,fullni[5].y,fullni[5].z};
    rij[0]=0; rij[1]=0; rij[2]=1;
    E+= bulk_potentialGPU(ni,nj,params,rij,1);
  }
  if (fullni[6].pt){
    float nj[3]={fullni[6].x,fullni[6].y,fullni[6].z};
    rij[0]=0; rij[1]=0; rij[2]=-1;
    E+= bulk_potentialGPU(ni,nj,params,rij,1);
  }

  return E;
}
__device__ float Bulk_latice_Potential_GPU(const nni fullni[], Parameters *params){
  float E=0;
  
  float ni[3]={fullni[0].x,fullni[0].y,fullni[0].z};
  
  E=newman_neighbours_GPU(fullni,params);
  if (params->elecA!=0)  {
    E+=Electric_Potential_GPU(ni,params)  ;
  }
  return E;
}

__device__ float Slab_latice_Potential_GPU(const nni fullni[], Parameters *params){
  float E=0;
  float ni[3]={fullni[0].x,fullni[0].y,fullni[0].z};
  
  E=newman_neighbours_GPU(fullni,params);
  
  float s[3]={fullni[7].x,fullni[7].y,fullni[7].z};
  if (fullni[0].pt>1) {
    E+= d_W[fullni[0].pt-2]*d_Surface_Potential_GPU[fullni[0].pt-2](ni,s,params,s);}
  if (params->elecA!=0)  {
    E+=Electric_Potential_GPU(ni,params)  ;
  }
  
  return E;
}

__device__ float Custom_latice_Potential_GPU(const nni fullni[], Parameters *params){
  float E=0;
  float ni[3]={fullni[0].x,fullni[0].y,fullni[0].z};
  
  E=newman_neighbours_GPU(fullni,params);
  float s[3]={fullni[7].x,fullni[7].y,fullni[7].z};

  if (fullni[0].pt>1) 
    E+= d_W[fullni[0].pt-2]*d_Surface_Potential_GPU[fullni[0].pt-2](ni,s,params,s);

  if (params->elecA!=0) 
    E+=Electric_Potential_GPU(ni,params);
      
  return E;
}
 
__device__ float Sphere_latice_Potential_GPU(const nni fullni[], Parameters *params){
  float E=0;
  float ni[3]={fullni[0].x,fullni[0].y,fullni[0].z};
  
  E=newman_neighbours_GPU(fullni,params);
  float s[3]={fullni[7].x,fullni[7].y,fullni[7].z};
  if (fullni[0].pt>1) 
    E+= d_W[fullni[0].pt-2]*d_Surface_Potential_GPU[fullni[0].pt-2](ni,s,params,s);

  if (params->elecA!=0) 
    E+=Electric_Potential_GPU(ni,params);
  
  return E;
}

__device__ float RP_Surface_Potential_GPU(float ni[3], float s[3], Parameters *params, float rij[3]){
  float nij=ni[0]*s[0]+ni[1]*s[1]+ni[2]*s[2];
  return -nij*nij;
}

__device__ float FG_Surface_Potential_GPU(float ni[3], float s[3], Parameters *params, float rij[3]){
  
  float nij=ni[0]*s[0]+ni[1]*s[1]+ni[2]*s[2];
  return +nij*nij;
}
__device__ float Strong_Anchoring_GPU(float ni[3], float s[3], Parameters *params, float rij[3]){

  return 0;
}

__device__ float RP_Surface_Potential_GHRL_GPU(float ni[3], float s[3], Parameters *params, float rij[3]){
  
  const float el=params->ghrl_lambda;
  const float em=params->ghrl_mu;
  const float en=params->ghrl_nu;
  const float er=params->ghrl_rho;
  const float es=params->ghrl_sigma;
  float v15 = 1.5;
  float v05 = 0.5;
  
  float ai=ni[0]*rij[0]+ni[1]*rij[1]+ni[2]*rij[2];
  float aj=s[0]*rij[0]+s[1]*rij[1]+s[2]*rij[2];
  float nij=ni[0]*s[0]+ni[1]*s[1]+ni[2]*s[2];
  float pij=v15*nij*nij-v05;
  float cross=(ni[2]*s[1]-ni[1]*s[2])*rij[0]
              +(ni[0]*s[2]-ni[2]*s[0])*rij[1]
              +(ni[1]*s[0]-ni[0]*s[1])*rij[2];
    float E1 =((v15*ai*ai)+(v15*aj*aj)-1);
  //~ if (threadIdx.x==0)printf("%0.3f %0.3f %0.3f %0.3f %0.3f \n", ai, aj, cross, nij,(ai*aj*nij));
  return((E1*(er*pij+el)+em*(ai*aj*nij)-(1/9))+en*pij+es*(nij>0?1:-1)*cross);
}
__device__ float FG_Surface_Potential_GHRL_GPU(float ni[3], float s[3], Parameters *params, float rij[3]){
  
  const float el=params->ghrl_lambda;
  const float em=params->ghrl_mu;
  const float en=params->ghrl_nu;
  const float er=params->ghrl_rho;
  const float es=params->ghrl_sigma;
  float v15 = 1.5;
  float v05 = 0.5;
  
  float ai=ni[0]*rij[0]+ni[1]*rij[1]+ni[2]*rij[2];
  float aj=s[0]*rij[0]+s[1]*rij[1]+s[2]*rij[2];
  float nij=ni[0]*s[0]+ni[1]*s[1]+ni[2]*s[2];
  float pij=v15*nij*nij-v05;
  float cross=(ni[2]*s[1]-ni[1]*s[2])*rij[0]
              +(ni[0]*s[2]-ni[2]*s[0])*rij[1]
              +(ni[1]*s[0]-ni[0]*s[1])*rij[2];
    float E1 =((v15*ai*ai)+(v15*aj*aj)-1);
  //~ if (threadIdx.x==0)printf("%0.3f %0.3f %0.3f %0.3f %0.3f \n", ai, aj, cross, nij,(ai*aj*nij));
  return -((E1*(er*pij+el)+em*(ai*aj*nij)-(1/9))+en*pij+es*(nij>0?1:-1)*cross);
}
__device__ float Strong_Anchoring_GHRL_GPU(float ni[3], float s[3], Parameters *params, float rij[3]){

  return 0;
}
