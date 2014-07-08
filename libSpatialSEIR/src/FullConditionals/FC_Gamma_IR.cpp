#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_FC_Gamma_IR.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    double dbeta(double x, double a, double b)
    {
        double out = (a-1)*std::log(x) + 
            (b-1)*std::log(1-x) + 
            (lgamma(a+b)) - 
            ((lgamma(a)) + (lgamma(b)));
        return(out);
    }



    FC_Gamma_IR::FC_Gamma_IR(ModelContext *_context,
                     CompartmentalModelMatrix *_R_star,
                     CompartmentalModelMatrix *_I,
                     InitData *_A0,
                     double *_p_ir,
                     double *_gamma_ir,
                     double _priorAlpha,
                     double _priorBeta,
                     int _useOCL,
                     double _sliceWidth)
    {

        context = new ModelContext*;
        R_star = new CompartmentalModelMatrix*;
        I = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_ir = new double*;
        gamma_ir = new double*;
        priorAlpha = new double;
        priorBeta = new double;
        value = new long double;
        samples = new int;
        accepted = new int; 
        useOCL = new int;
        sliceWidth = new double;
        *samples = 0;
        *accepted = 0;
        *useOCL = _useOCL;
        *context = _context;
        *R_star = _R_star;
        *I = _I;
        *A0 = _A0;
        *p_ir = _p_ir;
        *gamma_ir = _gamma_ir;
        *priorAlpha = _priorAlpha;
        *priorBeta = _priorBeta;
        *value = -1.0;
        *sliceWidth = _sliceWidth;

    }

    FC_Gamma_IR::~FC_Gamma_IR()
    {
        delete sliceWidth;
        delete gamma_ir;
        delete R_star;
        delete I;
        delete A0;
        delete p_ir;
        delete value;
        delete priorAlpha;
        delete priorBeta;
        delete context;
        delete samples;
        delete accepted;
        delete useOCL;
    }

    int FC_Gamma_IR::evalCPU()
    {
        *value = 0.0;
        int i, j, compIdx;
        int nLoc = *((*I) -> ncol);
        int nTpts = *((*I) -> nrow);
        for (i = 0; i < nLoc; i++)    
        {
            compIdx = i*nTpts;
            for (j = 0; j < nTpts; j++)     
            {
                *value += (*context) -> random -> dbinom(((*R_star) -> data)[compIdx],
                                                         ((*I) -> data)[compIdx],
                                                         (*p_ir)[j]);  
                compIdx++;
            }
        } 
        *value += (*context) -> random -> dgamma(**gamma_ir, *priorAlpha, 1/(*priorBeta));

        // Catch invalid values, nans etc. 
        if (!std::isfinite(*value))
        {
            *value = -INFINITY;
        }
        return(0);
    }

    int FC_Gamma_IR::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

    int FC_Gamma_IR::calculateRelevantCompartments()
    {
        (*context) -> calculateP_IR_CPU();
        return(0);
    }
    int FC_Gamma_IR::calculateRelevantCompartments_OCL()
    {
        return(calculateRelevantCompartments());
    }

    void FC_Gamma_IR::sample(int verbose)
    {
        if (verbose){std::cout << "Sampling P_IR\n";}
        if (*useOCL){sampleOCL(); return;}
        sampleCPU();
    }


    int FC_Gamma_IR::sampleCPU()
    {
        sampleEntireDouble_CPU(*context, *gamma_ir, 1, *sliceWidth); 
        return(0);
    }

    int FC_Gamma_IR::sampleOCL()
    {
        //NOT IMPLEMENTED
        return(sampleCPU());
    }

    long double FC_Gamma_IR::getValue()
    {
        return(*(this -> value));
    }
    void FC_Gamma_IR::setValue(long double val)
    {
        *(this -> value) = val;
    }
}
