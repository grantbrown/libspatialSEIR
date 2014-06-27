#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_FC_Rho.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;


    FC_Rho::FC_Rho(ModelContext *_context,
                   CompartmentalModelMatrix *_E_star,  
                   CompartmentalModelMatrix *_S,
                   InitData *_A0,
                   CovariateMatrix *_X,
                   double *_p_se,
                   double *_beta,
                   double *_rho,
                   double _sliceWidth)
    {
        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        X = new CovariateMatrix*;
        p_se = new double*;
        beta = new double*;
        rho = new double*;
        sliceWidth = new double;
        value = new long double;
        samples = new int;
        accepted = new int; 
        *samples = 0;
        *accepted = 0;


        *context = _context;
        *E_star = _E_star;
        *S = _S;
        *A0 = _A0;
        *X = _X;
        *p_se = _p_se;
        *beta = _beta;
        *rho = _rho;
        *sliceWidth = _sliceWidth;
        *value = -1.0;
    }
    FC_Rho::~FC_Rho()
    {
        delete E_star;
        delete S;
        delete A0;
        delete X;
        delete p_se;
        delete beta;
        delete rho;
        delete value;
        delete sliceWidth;
        delete context;
        delete samples;
        delete accepted;

    }

    int FC_Rho::evalCPU()
    {
        *value = 0.0;
        int i, j, Es, compIdx;
        double pse;
        int nLoc = *((*S) -> ncol);
        int nTpts = *((*S) -> nrow);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;

        for (i = 0; i < nLoc; i++)    
        {
            compIdx = i*nTpts;
            for (j = 0; j < nTpts; j++)     
            {
                // todo: use dbinom
                Es = ((*E_star) -> data)[compIdx];
                pse = (*p_se)[compIdx];
                term1 += std::log(pse)*Es; 
                term2 += std::log(1-pse)*(((*S) -> data)[compIdx] - Es);
                compIdx++;
            }
        } 
        term3 += (**rho > 0 && **rho < 1 ? 0 : -INFINITY); // Generalize to allow informative priors. 
                                                        // Prior specification in this area needs work. 
        *value = term1 + term2 + term3;
        // Catch invalid values, nans etc. 
        if (!std::isfinite(*value))
        {
            *value = -INFINITY;
        }

        return(0);
    }

    int FC_Rho::evalOCL()
    {
        //NOT IMPLEMENTED
        return(evalCPU());
    }
    int FC_Rho::calculateRelevantCompartments()
    {
       (*context) -> calculateP_SE_CPU();
       return(0); 
    }
    int FC_Rho::calculateRelevantCompartments_OCL()
    {
       (*context) -> calculateP_SE_OCL();
       return(0); 
    }

    int FC_Rho::sampleCPU()
    {
        sampleDoubleMetropolis(*context, *rho, 1, *sliceWidth); 
        return(0);
    }
    int FC_Rho::sampleOCL()
    {
        sampleEntireDouble_OCL(*context, *rho, 1, *sliceWidth); 
        return(0);
    }

    long double FC_Rho::getValue()
    {
        return(*(this -> value));
    }
    void FC_Rho::setValue(long double val)
    {
        *(this -> value) = val;
    }
}
