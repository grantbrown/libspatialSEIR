#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_FC_P_EI.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>


double dbeta(double x, double a, double b)
{
    double out = (a-1)*std::log(x) + 
        (b-1)*std::log(1-x) + 
        (lgamma(a+b)) - 
        ((lgamma(a)) + (lgamma(b)));
    return(out);
}

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;



    FC_P_EI::FC_P_EI(ModelContext *_context,
                     CompartmentalModelMatrix *_I_star,
                     CompartmentalModelMatrix *_E,
                     InitData *_A0,
                     double *_p_ei,
                     double _priorAlpha,
                     double _priorBeta,
                     int _useOCL)
    {

        context = new ModelContext*;
        I_star = new CompartmentalModelMatrix*;
        E = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_ei = new double*;
        priorAlpha = new double;
        priorBeta = new double;
        value = new long double;
        samples = new int;
        accepted = new int; 
        useOCL = new int;
        *samples = 0;
        *accepted = 0;
        *useOCL = _useOCL;
        *context = _context;
        *I_star = _I_star;
        *E = _E;
        *A0 = _A0;
        *p_ei = _p_ei;
        *priorAlpha = _priorAlpha + 1;
        *priorBeta = _priorBeta + 1;
        *value = -1.0;

    }
    FC_P_EI::~FC_P_EI()
    {
        delete I_star;
        delete E;
        delete A0;
        delete p_ei;
        delete value;
        delete priorAlpha;
        delete priorBeta;
        delete context;
        delete samples;
        delete accepted;
        delete useOCL;
    }

    int FC_P_EI::evalCPU()
    { 
        *value = 0.0;
        int i_star_sum = (*I_star) -> marginSum(3,-1);
        int e_sum = (*E) -> marginSum(3,-1);
        *value = dbeta(**p_ei, *priorAlpha + i_star_sum, *priorBeta - i_star_sum + e_sum); 
        return 0;
    }

    int FC_P_EI::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_EI::calculateRelevantCompartments()
    {
        // Not used, Do nothing
        return(0);
    }
    int FC_P_EI::calculateRelevantCompartments_OCL()
    {
        // Not used, Do nothing
        return(0);
    }

    void FC_P_EI::sample(int verbose)
    {
        if (verbose){std::cout << "Sampling P_EI\n";}
        if (*useOCL){sampleOCL(); return;}
        sampleCPU();
    }

    int FC_P_EI::sampleCPU()
    {
        double a, b;
        a = ((*I_star) -> marginSum(3, -1));
        b = ((*E) -> marginSum(3, -1)) - a;
        //std::cout << "(a,b): (" << a << "," << b <<")\n";
        //std::cout << "(a,b): (" << a + *priorAlpha << "," << b+*priorBeta <<")\n";
        (**p_ei) = ((*context) -> random -> beta(a+*priorAlpha, b+*priorBeta));
        return(0);
    }
    int FC_P_EI::sampleOCL()
    {
        //NOT IMPLEMENTED
        return(sampleCPU());
    }

    long double FC_P_EI::getValue()
    {
        return(*(this -> value));
    }
    void FC_P_EI::setValue(long double val)
    {
        *(this -> value) = val;
    }
}
