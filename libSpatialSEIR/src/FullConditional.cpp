#ifndef SPATIALSEIR_INCLUDEFILES
#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#endif
#include<cblas.h>
#include<cmath>
#ifndef FULL_CONDITIONAL_INC
#define FULL_CONDITIONAL_INC
#include <FullConditional.hpp>
#endif

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    /*
     * Helper functions
     */
    int matMult(double* output, double * A, double * B, int Arow, int Acol, int Brow, int Bcol, bool TransA = false, bool TransB = false)
    {
        // Use BLAS to matrix multiply, assume column major, non transposing.
        // Double check this code when I have internet access. 
        cblas_dgemm(CblasColMajor,
                    TransA ? CblasTrans : CblasNoTrans,
                    TransB ? CblasTrans : CblasNoTrans,
                    Arow, Bcol, Brow,
                    1.0, 
                    A, Arow, 
                    B, Brow, 
                    0.0, output, Brow);
        return 0; 
    }

    double dbeta(double x, double a, double b)
    {
        // Not Implemented
        return(0.0);
    }


    /*
     *
     * Implement the data container class InitData
     *
     */    
 
 
    InitData::InitData(int *_S0, 
                       int *_E0,
                       int *_I0,
                       int *_R0,
                       int *_S_star0,
                       int *_E_star0, 
                       int *_I_star0, 
                       int *_R_star0,
                       int *nLoc)
    {
        this -> populate(*&_S0, *&_E0, *&_I0, *&_R0, 
                *&_S_star0, *&_E_star0, *&_I_star0, *&_R_star0,
                *&nLoc);
    }

    InitData::InitData()
    {
        // Do nothing
    }

    void InitData::populate(int *_S0, 
                       int *_E0,
                       int *_I0,
                       int *_R0,
                       int *_S_star0,
                       int *_E_star0, 
                       int *_I_star0, 
                       int *_R_star0,
                       int *nLoc
                       )
    {
        S0 = new int[*nLoc];
        E0 = new int[*nLoc];
        I0 = new int[*nLoc];
        R0 = new int[*nLoc];
        S_star0 = new int[*nLoc]; 
        E_star0 = new int[*nLoc];
        I_star0 = new int[*nLoc];
        R_star0 = new int[*nLoc];
        numLocations = new int;
        *numLocations = *nLoc;
        int i;
        for (i = 0; i < *nLoc; i++)
        {
            S0[i] = _S0[i]; 
            E0[i] = _E0[i];
            I0[i] = _I0[i];
            R0[i] = _R0[i];
            S_star0[i] = _S_star0[i];
            E_star0[i] = _E_star0[i];
            I_star0[i] = _I_star0[i];
            R_star0[i] = _R_star0[i];
        }
    }

    InitData::~InitData()
    {
        delete S0;
        delete E0;
        delete I0;
        delete R0;
        delete S_star0;
        delete E_star0;
        delete I_star0;
        delete R_star0;
        delete numLocations;
    }

    FullConditional::FullConditional()
    {
        //This class should not be instantiated directly. 
        throw(-1);
    }
    FullConditional::~FullConditional()
    {
        throw(-1);
    }
    int FullConditional::evalCPU()
    {
        throw(-1);
    }
    int FullConditional::evalOCL()
    {
        throw(-1);
    }
    int FullConditional::sampleCPU()
    {
        throw(-1);
    }
    int FullConditional::sampleOCL()
    {
        throw(-1);
    }
     
    /*
     *
     * Implement the full conditional distribution for S_star
     *
     */    

    FC_S_Star::FC_S_Star(ModelContext * _context,
                         CompartmentalModelMatrix *_S_star, 
                         CompartmentalModelMatrix *_S, 
                         CompartmentalModelMatrix *_R, 
                         InitData *_A0,
                         CovariateMatrix *_X, 
                         double *_p_se,
                         double *_p_rs,
                         double *_beta,
                         double *_rho)
    {
       context = new ModelContext*;
       S_star = new CompartmentalModelMatrix*;
       S = new CompartmentalModelMatrix*;
       R = new CompartmentalModelMatrix*;
       A0 = new InitData*;
       X = new CovariateMatrix*;
       p_se = new double*;
       p_rs = new double*;
       beta = new double*;
       rho = new double*;
       value = new double;
       *context = _context;
       *S_star = _S_star;
       *S = _S;
       *R = _R;
       *A0 = _A0;
       *X = _X;
       *p_se = _p_se;
       *p_rs = _p_rs;
       *beta = _beta;
       *rho = _rho;
       *value = -1.0;
    }    
    FC_S_Star::~FC_S_Star()
    {
        delete S_star;
        delete S;
        delete R;
        delete A0;
        delete X;
        delete p_se;
        delete p_rs;
        delete beta;
        delete rho;
        delete value;
        delete context;
    }

    // Evaluate the S_star FC at the current values provided by the context.
    int FC_S_Star::evalCPU()
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        for (j = 0; j < nTpts; j++)     
        {
            for (i = 0; i < nLoc; i++)    
            {
                compIdx = i + j*nLoc;
                tmp = ((*S_star) -> data)[compIdx];
                if (tmp < 0)
                {
                    *value = -std::log(0.0);
                    return(-1);
                }
                term1 += std::log((*p_rs)[j])*tmp; 
                term2 += std::log(1-(*p_rs)[j])*(((*R) -> data)[compIdx] - tmp);
                term3 += std::log(1-(*p_se)[compIdx])*(((*S) -> data)[compIdx]) ;
            }
        } 
        *value = term1 + term2 + term3;
        return(0);
    }
    int FC_S_Star::evalOCL()
    {
        //NOT IMPLEMENTED
        return-1;
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

    
    FC_E_Star::FC_E_Star(ModelContext *_context,
                         CompartmentalModelMatrix *_E_star,
                         CompartmentalModelMatrix *_E,  
                         CompartmentalModelMatrix *_S,
                         CovariateMatrix *_X,
                         InitData *_A0,
                         double *_p_se,
                         double *_p_ei,
                         double *_rho,
                         double *_beta) 
    {

        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        E = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        X = new CovariateMatrix*;
        A0 = new InitData*;
        p_se = new double*;
        p_ei = new double*;
        rho = new double*;
        beta = new double*;
        value = new double;
        
       *context = _context;
        *E_star = _E_star;
        *E = _E;
        *S = _S;
        *X = _X;
        *A0 = _A0;
        *p_se = _p_se;
        *p_ei = _p_ei;
        *rho = _rho;
        *beta = _beta;
        *value = -1.0;
    }

    FC_E_Star::~FC_E_Star()
    {
        delete E_star;
        delete E;
        delete S;
        delete X;
        delete A0;
        delete p_se;
        delete p_ei;
        delete rho;
        delete beta;
        delete value;
        delete context;
    }

    int FC_E_Star::evalCPU()
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        for (j = 0; j < nTpts; j++)     
        {
            for (i = 0; i < nLoc; i++)    
            {
                compIdx = i + j*nLoc;
                tmp = ((*E_star) -> data)[compIdx];
                if (tmp < 0)
                {
                    *value = -std::log(0.0);
                    return(-1);
                }
                term1 += std::log((*p_se)[compIdx])*tmp; 
                term2 += std::log(1-(*p_se)[compIdx])*(((*S) -> data)[compIdx] - tmp);
                term3 += std::log(1-(**p_ei))*(((*E) -> data)[compIdx]) ;
            }
        } 
        *value = term1 + term2 + term3;
        return(0);
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
    FC_R_Star::FC_R_Star(ModelContext *_context,
                         CompartmentalModelMatrix *_R_star,
                         CompartmentalModelMatrix *_R,
                         CompartmentalModelMatrix *_I,
                         InitData *_A0,
                         double *_p_rs,
                         double *_p_ir)
    {

        context = new ModelContext*;
        R_star = new CompartmentalModelMatrix*;
        R = new CompartmentalModelMatrix*;
        I = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_rs = new double*;
        p_ir = new double*;
        value = new double;

       *context = _context;
        *R_star = _R_star;
        *R = _R;
        *I = _I;
        *A0 = _A0;
        *p_rs = _p_rs;
        *p_ir = _p_ir;
        *value = -1.0;
    }
    FC_R_Star::~FC_R_Star()
    {
        delete R_star;
        delete R;
        delete I;
        delete A0;
        delete p_rs;
        delete p_ir;
        delete value;
        delete context;
    }

    int FC_R_Star::evalCPU()
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*I) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        for (j = 0; j < nTpts; j++)     
        {
            for (i = 0; i < nLoc; i++)    
            {
                compIdx = i + j*nLoc;
                tmp = ((*R_star) -> data)[compIdx];
                if (tmp < 0)
                {
                    *value = -std::log(0.0);
                    return(-1);
                }
                term1 += std::log((**p_ir))*tmp; 
                term2 += std::log(1-(**p_ir))*(((*I) -> data)[compIdx] - tmp);
                term3 += std::log(1-((*p_rs)[j]))*(((*R) -> data)[compIdx]) ;
            }
        } 
        *value = term1 + term2 + term3;
        return(0);
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


    FC_Beta::FC_Beta(ModelContext *_context,
                     CompartmentalModelMatrix *_E_star, 
                     CompartmentalModelMatrix *_S, 
                     InitData *_A0,
                     CovariateMatrix *_X,
                     double *_p_se, 
                     double *_beta, 
                     double *_rho)
    {

        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        X = new CovariateMatrix*;
        p_se = new double*;
        beta = new double*;
        rho = new double*;
        value = new double;

        *context = _context;
        *E_star = _E_star;
        *S = _S;
        *A0 = _A0;
        *X = _X;
        *p_se = _p_se;
        *beta = _beta;
        *rho = _rho;
        *value = -1.0;
    }

    FC_Beta::~FC_Beta()
    {
        delete E_star;
        delete S;
        delete A0;
        delete X;
        delete p_se;
        delete beta;
        delete rho;
        delete value;
        delete context;
    }
    
    int FC_Beta::evalCPU()
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        for (j = 0; j < nTpts; j++)     
        {
            for (i = 0; i < nLoc; i++)    
            {
                compIdx = i + j*nLoc;
                tmp = ((*E_star) -> data)[compIdx];
                if (tmp < 0)
                {
                    *value = -std::log(0.0);
                    return(-1);
                }
                term1 += std::log((*p_se)[compIdx])*tmp; 
                term2 += std::log(1-(*p_se)[compIdx])*(((*S) -> data)[compIdx] - tmp);
            }
        } 
        for (i = 0; i < (*((*X) -> ncol_x) + *((*X) -> ncol_z)); i++)
        {
            term3 += pow((*beta)[i],2)/10; // Generalize to allow different prior precisions. 
        }
        *value = term1 + term2 + term3;
        return(0);
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
    FC_P_RS::FC_P_RS(ModelContext *_context,
                     CompartmentalModelMatrix *_S_star, 
                     CompartmentalModelMatrix *_R,
                     InitData *_A0,
                     double *_p_rs)
    {

        context = new ModelContext*;
        S_star = new CompartmentalModelMatrix*;
        R = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_rs = new double*;
        value = new double;

        *context = _context;
        *S_star = _S_star;
        *R = _R;
        *A0 = _A0;
        *p_rs = _p_rs;
        *value = -1.0;
    }
    FC_P_RS::~FC_P_RS()
    {
        delete S_star;
        delete R;
        delete A0;
        delete p_rs;
        delete value;
        delete context;
    }
    int FC_P_RS::evalCPU()
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*R) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        int* s_star_i_sum = new int[nTpts];
        int* r_star_i_diff = new int[nTpts];
        for (j =0; j<nTpts; j++)
        {
            s_star_i_sum[j]=0.0;
            r_star_i_diff[j]=0.0;
        }


        for (j = 0; j < nTpts; j++)     
        {
            for (i = 0; i < nLoc; i++)    
            {
               compIdx = i + j*nLoc;
                tmp = ((*S_star) -> data)[compIdx];
                if (tmp < 0)
                {
                    *value = -std::log(0.0);
                    return(-1);
                }
                s_star_i_sum[j] += tmp; 
                r_star_i_diff[j] += (((*R) -> data)[compIdx] - tmp);
            }
        } 
        *value = 0;
        for (j = 0; j < nTpts; j++)
        {
            *value += dbeta((*p_rs)[j],1.5 + s_star_i_sum[j], 1.5 + r_star_i_diff[j]);
        }
       
        delete[] s_star_i_sum;
        delete[] r_star_i_diff;
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


    FC_Rho::FC_Rho(ModelContext *_context,
                   CompartmentalModelMatrix *_S_star,
                   CompartmentalModelMatrix *_E_star,
                   InitData *_A0,
                   CovariateMatrix *_X,
                   double *_p_se,
                   double *_beta,
                   double *_rho)
    {
        context = new ModelContext*;
        S_star = new CompartmentalModelMatrix*;
        E_star = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        X = new CovariateMatrix*;
        p_se = new double*;
        beta = new double*;
        rho = new double*;
        value = new double;

        *context = _context;
        *S_star = _S_star;
        *E_star = _E_star;
        *A0 = _A0;
        *X = _X;
        *p_se = _p_se;
        *beta = _beta;
        *rho = _rho;
        *value = -1.0;
    }
    FC_Rho::~FC_Rho()
    {
        delete S_star;
        delete E_star;
        delete A0;
        delete X;
        delete p_se;
        delete beta;
        delete rho;
        delete value;
        delete context;
    }
    int FC_Rho::evalCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_Rho::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_Rho::sampleCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_Rho::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }




    FC_P_EI::FC_P_EI(ModelContext *_context,
                     CompartmentalModelMatrix *_I_star,
                     CompartmentalModelMatrix *_E_star,
                     InitData *_A0,
                     double *_p_ei)
    {

        context = new ModelContext*;
        I_star = new CompartmentalModelMatrix*;
        E_star = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_ei = new double*;
        value = new double;

        *context = _context;
        *I_star = _I_star;
        *E_star = _E_star;
        *A0 = _A0;
        *p_ei = _p_ei;
        *value = -1.0;

    }
    FC_P_EI::~FC_P_EI()
    {
        delete I_star;
        delete E_star;
        delete A0;
        delete p_ei;
        delete value;
        delete context;
    }
    int FC_P_EI::evalCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_EI::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_EI::sampleCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_EI::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }



    FC_P_IR::FC_P_IR(ModelContext *_context,
                     CompartmentalModelMatrix *_I_star,
                     CompartmentalModelMatrix *_R_star,
                     InitData *_A0,
                     double *_p_ir)
    {

        context = new ModelContext*;
        I_star = new CompartmentalModelMatrix*;
        R_star = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_ir = new double*;
        value = new double;

        *context = _context;
        *I_star = _I_star;
        *R_star = _R_star;
        *A0 = _A0;
        *p_ir = _p_ir;
        *value = -1.0;

    }
    FC_P_IR::~FC_P_IR()
    {
        delete I_star;
        delete R_star;
        delete A0;
        delete p_ir;
        delete value;
        delete context;
    }
    int FC_P_IR::evalCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_IR::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_IR::sampleCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_IR::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
}


