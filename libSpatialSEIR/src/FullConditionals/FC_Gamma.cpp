#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_FC_Gamma.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    // Depricated: to be removed
    FC_Gamma::FC_Gamma(ModelContext *_context,
                   CompartmentalModelMatrix *_E_star,  
                   CompartmentalModelMatrix *_S,
                   InitData *_A0,
                   CovariateMatrix *_X,
                   double *_p_se,
                   double *_beta,
                   double *_gamma,
                   double *_priorAlpha,
                   double *_priorBeta,
                   double _sliceWidth,
                   int _useOCL)
    {
        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        X = new CovariateMatrix*;
        p_se = new double*;
        beta = new double*;
        gamma = new double*;
        priorAlpha = new double;
        priorBeta = new double;
        sliceWidth = new double;
        value = new long double;
        samples = new int;
        accepted = new int; 
        *samples = 0;
        *accepted = 0;
        useOCL = new int;



        *context = _context;
        *E_star = _E_star;
        *S = _S;
        *A0 = _A0;
        *X = _X;
        *p_se = _p_se;
        *beta = _beta;
        *gamma = _gamma;
        *priorAlpha = *_priorAlpha;
        *priorBeta = *_priorBeta;
        *sliceWidth = _sliceWidth;
        *value = -1.0;
        *useOCL = _useOCL;
    }
    FC_Gamma::~FC_Gamma()
    {
        delete E_star;
        delete S;
        delete A0;
        delete X;
        delete p_se;
        delete beta;
        delete gamma;
        delete priorAlpha;
        delete priorBeta;
        delete sliceWidth;
        delete value;
        delete context;
        delete samples;
        delete accepted;
        delete useOCL;

    }

    int FC_Gamma::evalCPU()
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
                Es = ((*E_star) -> data)[compIdx];
                pse = (*p_se)[compIdx];
                term1 += std::log(pse)*Es; 
                term2 += std::log(1-pse)*(((*S) -> data)[compIdx] - Es);
                compIdx++;
            }
        } 
        for (j = 0; j < nTpts; j++)
        {
            term3 += ((*priorAlpha-1)*std::log((*gamma)[j]) - ((*gamma)[j])/(*priorBeta)); 
        }
        *value = term1 + term2 + term3;
        // Catch invalid values, nans etc. 
        if (!std::isfinite(*value))
        {
            *value = -INFINITY;
        }

        return(0);
    }

    int FC_Gamma::evalOCL()
    {
        //NOT IMPLEMENTED
        return(evalCPU());
    }
    int FC_Gamma::calculateRelevantCompartments()
    {
       (*context) -> calculateP_SE_CPU();
       return(0); 
    }
    int FC_Gamma::calculateRelevantCompartments_OCL()
    {
       (*context) -> calculateP_SE_OCL();
       return(0); 
    }

    void FC_Gamma::sample(int verbose)
    {
        if (verbose){std::cout << "Sampling Gamma\n";}
        if (*useOCL){sampleOCL(); return;}
        sampleCPU();
    }

    int FC_Gamma::sampleCPU()
    {
        sampleEntireDouble_CPU(*context, *gamma, *((*A0) -> numLocations), sliceWidth); 
        return(0);
    }
    int FC_Gamma::sampleOCL()
    {
        sampleEntireDouble_OCL(*context, *gamma, *((*A0) -> numLocations), sliceWidth); 
        return(0);
    }

    long double FC_Gamma::getValue()
    {
        return(*(this -> value));
    }
    void FC_Gamma::setValue(long double val)
    {
        *(this -> value) = val;
    }
}
