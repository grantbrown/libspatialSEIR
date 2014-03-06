
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

    /*
     *
     * Implement the data container class InitData
     *
     */    
 
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
     
    /*
     *
     * Implement the full conditional distribution for S_star
     *
     */    

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


    /*
     *
     * Implement the full conditional distribution for E_Star
     *
     */    

    
    FC_E_Star::FC_E_Star(CompartmentalModelMatrix *_S_star,
                         CompartmentalModelMatrix *_E_star,
                         CompartmentalModelMatrix *_I_star,
                         CovariateMatrix *_X,
                         InitData *_A0,
                         double *_p_se,
                         double *_p_rs,
                         double *_rho,
                         double *_beta) 
    {
        CompartmentalModelMatrix* S_star = _S_star;
        CompartmentalModelMatrix* E_star = _E_star;
        CompartmentalModelMatrix* I_star = _I_star;
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

    /*
     *
     * Implement the full conditional distribution for R_Star
     *
     */    
    FC_R_Star::FC_R_Star(CompartmentalModelMatrix *_R_star,
                         CompartmentalModelMatrix *_S_star,
                         CompartmentalModelMatrix *_I_star,
                         InitData *_A0,
                         double *_p_rs,
                         double *_p_ir)
    {
        CompartmentalModelMatrix* R_star = _R_star;
        CompartmentalModelMatrix* S_star = _S_star;
        CompartmentalModelMatrix* I_star = _I_star;
        InitData* A0 = _A0;
        double* p_rs = _p_rs;
        double* p_ir = _p_ir;
    }

    int FC_R_Star::evalCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_R_Star::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_R_Star::sampleCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_R_Star::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

     /*
     *
     * Implement the full conditional distribution for the regression
     * parameters: beta
     */


    FC_Beta::FC_Beta(CompartmentalModelMatrix *_E_star, 
                     CompartmentalModelMatrix *_S_star, 
                     InitData *_A0,
                     CovariateMatrix *_X,
                     double *_p_se, 
                     double *_beta, 
                     double *_rho)
    {
        CompartmentalModelMatrix* E_star = _E_star;
        CompartmentalModelMatrix* S_star = _S_star;
        InitData* A0 = _A0;
        CovariateMatrix* X = _X;
        double* p_se = _p_se;
        double* beta = _beta;
        double* rho = _rho;
    }
    int FC_Beta::evalCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_Beta::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_Beta::sampleCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_Beta::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    /*
     *
     * Implement the full conditional for the R->S transition 
     * probabilities. 
     *
     */

    FC_P_RS::FC_P_RS(CompartmentalModelMatrix *_S_star, 
                     CompartmentalModelMatrix *_R_star,
                     InitData *_A0,
                     double *_p_rs)
    {
        CompartmentalModelMatrix* S_star = _S_star;
        CompartmentalModelMatrix* R_star = _R_star;
        InitData* A0 = _A0;
        double* p_rs = _p_rs;
    }
    int FC_P_RS::evalCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_RS::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_RS::sampleCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_RS::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }


    FC_Rho::FC_Rho(CompartmentalModelMatrix *_S_star,
                   CompartmentalModelMatrix *_E_star,
                   InitData *_A0,
                   CovariateMatrix *_X,
                   double *_p_se,
                   double *_beta,
                   double *_rho)
    {
        CompartmentalModelMatrix* S_star = _S_star;
        CompartmentalModelMatrix* E_star = _E_star;
        InitData* A0 = _A0;
        CovariateMatrix* X = _X;
        double* p_se = _p_se;
        double* beta = _beta;
        double* rho = rho;
    }
}


