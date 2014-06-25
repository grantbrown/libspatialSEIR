#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<FullConditional.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;


    /*
     *
     * Implement the full conditional for the R->S transition 
     * probabilities. 
     *
     */
    FC_Beta_P_RS::FC_Beta_P_RS(ModelContext *_context,
                     CompartmentalModelMatrix *_S_star, 
                     CompartmentalModelMatrix *_R,
                     CovariateMatrix* _X,
                     InitData *_A0,
                     double *_p_rs,
                     double *_beta_p_rs,
                     double _tausq,
                     double _sliceWidth)
    {

        context = new ModelContext*;
        S_star = new CompartmentalModelMatrix*;
        R = new CompartmentalModelMatrix*;
        X = new CovariateMatrix*;
        A0 = new InitData*;
        p_rs = new double*;
        beta_p_rs = new double*;
        tausq = new double;
        sliceWidth = new double;
        value = new long double;
        samples = new int; *samples = 0;
        accepted = new int; *accepted = 0;

        *context = _context;
        *S_star = _S_star;
        *X = _X;
        *R = _R;
        *A0 = _A0;
        *p_rs = _p_rs;
        *beta_p_rs = _beta_p_rs;
        *tausq = _tausq;
        *sliceWidth = _sliceWidth;
        *value = -1.0;
    }
    FC_Beta_P_RS::~FC_Beta_P_RS()
    {
        delete S_star;
        delete R;
        delete X;
        delete beta_p_rs;
        delete tausq;
        delete A0;
        delete p_rs;
        delete value;
        delete sliceWidth;
        delete context;
        delete samples;
        delete accepted;

    }

    int FC_Beta_P_RS::evalCPU()
    {
        int j;
        long double a,b;
        int nbeta = *((*X) -> ncol_x);
        int nTpts = *((*R) -> nrow);
        double tmp;
        long double term1 = 0.0;
        double term2 = 0.0;


        for (j = 0; j < nTpts; j++)
        {
            tmp = (*p_rs)[j];
            if (tmp <= 0 || tmp >= 1)
            {
                *value = -INFINITY;
                return(-1);
            }
            a = ((*S_star)-> marginSum(1,j));
            b = ((*R) -> marginSum(1,j)); 
            term1 += std::log(tmp)*(a);
            term1 += std::log(1-tmp)*(b-a);
        }

        for (j = 0; j < nbeta; j++)
        {
            term2 -= ((*tausq)/2)*pow((*beta_p_rs)[j],2);
        }
        *value = term1 + term2;
        if (!std::isfinite(*value))
        {
            *value = -INFINITY;
        }
        return(0);
    }

    int FC_Beta_P_RS::evalOCL()
    {
        //NOT IMPLEMENTED
        return(evalCPU());
    }
    int FC_Beta_P_RS::calculateRelevantCompartments()
    {
         ((*context) -> calculateP_RS_CPU());      
         return(0);
    }
    int FC_Beta_P_RS::calculateRelevantCompartments_OCL()
    {
         ((*context) -> calculateP_RS_CPU());      
         return(0);
    }
    int FC_Beta_P_RS::sampleCPU()
    {
        int nbeta = *((*X) -> ncol_x);
        int mode = (*context) -> getSamplingMode();

        if (mode == 1)
        {
            sampleDoubleMetropolis(*context, *beta_p_rs, nbeta, *sliceWidth); 
        }
        else
        {
            sampleDoubleMetropolis(*context, *beta_p_rs, nbeta, *sliceWidth); 
        }
        return(0);
    }

    int FC_Beta_P_RS::sampleOCL()
    {
        int nbeta = *((*X) -> ncol_x);
        sampleEntireDouble_OCL(*context, *beta_p_rs, nbeta, *sliceWidth); 
        return(0);
    }

    long double FC_Beta_P_RS::getValue()
    {
        return(*(this -> value));
    }
    void FC_Beta_P_RS::setValue(long double val)
    {
        *(this -> value) = val;
    }
}
