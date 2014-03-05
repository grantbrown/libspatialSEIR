
#ifndef SPATIALSEIR_INCLUDEFILES
#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#endif

#include "FullConditional.hpp"



namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    //Implement data containder class InitData
    
    InitData::InitData(double *_S0, 
                       double *_E0,
                       double *_I0,
                       double *_R0,
                       double *_S_star0,
                       double *_E_star0, 
                       double *_I_star0, 
                       double *_R_star0)
    {
        double* S0 = _S0;
        double* E0 = _E0;
        double* I0 = _I0;
        double* R0 = _R0;
        double* S_star0 = _S_star0;
        double* E_star0 = _E_star0;
        double* I_star0 = _I_star0;
        double* R_star0 = _R_star0;
    }

    //Implement FC_S_Star
     
    FC_S_Star::FC_S_Star(CompartmentalModelMatrix *_S_star, 
                         CompartmentalModelMatrix *_E_star, 
                         CompartmentalModelMatrix *_R_star, 
                         InitData *_A0,
                         CovariateMatrix *_X, 
                         double *_p_se,
                         double *_p_rs,
                         double *_beta,
                         double *_rho)
    {
       CompartmentalModelMatrix* S_star = _S_star;
       CompartmentalModelMatrix* E_star = _E_star;
       CompartmentalModelMatrix* R_star = _R_star;
       InitData* A0 = _A0;
       CovariateMatrix*X = _X;
       double* p_se = _p_se;
       double* p_rs = _p_rs;
       double* beta = _beta;
       double* rho = _rho;
       double* value; *value = -1.0;
    }    

    int FC_S_Star::evalCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_S_Star::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_S_Star::sampleCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_S_Star::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

    //Implement FC_E_Star
    
    FC_E_Star::FC_E_Star(CompartmentalModelMatrix *_S_Star,
                         CompartmentalModelMatrix *_E_Star,
                         CompartmentalModelMatrix *_I_Star,
                         CovariateMatrix *_X,
                         InitData *_A0,
                         double *_p_se,
                         double *_p_rs,
                         double *_rho,
                         double *_beta) 
    {
        CompartmentalModelMatrix* S_Star = _S_Star;
        CompartmentalModelMatrix* E_Star = _E_Star;
        CompartmentalModelMatrix* I_Star = _I_Star;
        CovariateMatrix* X = _X;
        InitData* A0 = _A0;
        double* p_se = _p_se;
        double* p_rs = _p_rs;
        double* rho = _rho;
        double* beta = _beta;
        double* value; *value = -1.0;
    }

    int FC_E_Star::evalCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_E_Star::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_E_Star::sampleCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_E_Star::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }


}
