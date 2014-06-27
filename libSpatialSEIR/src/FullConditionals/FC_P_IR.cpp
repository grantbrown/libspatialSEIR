#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_FC_P_IR.hpp>
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



    FC_P_IR::FC_P_IR(ModelContext *_context,
                     CompartmentalModelMatrix *_R_star,
                     CompartmentalModelMatrix *_I,
                     InitData *_A0,
                     double *_p_ir,
                     double _priorAlpha,
                     double _priorBeta)
    {

        context = new ModelContext*;
        R_star = new CompartmentalModelMatrix*;
        I = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_ir = new double*;
        priorAlpha = new double;
        priorBeta = new double;
        value = new long double;
        samples = new int;
        accepted = new int; 
        *samples = 0;
        *accepted = 0;

        *context = _context;
        *R_star = _R_star;
        *I = _I;
        *A0 = _A0;
        *p_ir = _p_ir;
        *priorAlpha = _priorAlpha + 1;
        *priorBeta = _priorBeta + 1;
        *value = -1.0;

    }

    FC_P_IR::~FC_P_IR()
    {
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
    }

    int FC_P_IR::evalCPU()
    {
        *value = 0.0;
        int r_star_sum = (*R_star) -> marginSum(3,-1);
        int i_sum = (*I) -> marginSum(3,1);
        *value = dbeta(**p_ir, *priorAlpha + r_star_sum, *priorBeta - r_star_sum + i_sum); 
        return 0;
    }

    int FC_P_IR::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

    int FC_P_IR::calculateRelevantCompartments()
    {
        // not used, do nothing. 
        return(0);
    }
    int FC_P_IR::calculateRelevantCompartments_OCL()
    {
        // Not used, Do nothing
        return(0);
    }

    int FC_P_IR::sampleCPU()
    {
        double a,b;
        a = (*R_star) -> marginSum(3,-1);
        b = ((*I) -> marginSum(3,-1)) - a;
        (**p_ir) = ((*context)->random->beta(a+(*priorAlpha), b+(*priorBeta)));
        return(0);
    }

    int FC_P_IR::sampleOCL()
    {
        //NOT IMPLEMENTED
        return(sampleCPU());
    }

    long double FC_P_IR::getValue()
    {
        return(*(this -> value));
    }
    void FC_P_IR::setValue(long double val)
    {
        *(this -> value) = val;
    }
}
