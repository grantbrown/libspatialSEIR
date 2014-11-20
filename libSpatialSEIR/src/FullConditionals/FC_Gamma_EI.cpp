#include<math.h>
#include<cstring>
#include<vector>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FC_Gamma_EI.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>
#include<IOProvider.hpp>

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

    FC_Gamma_EI::FC_Gamma_EI(ModelContext *_context,
                     CompartmentalModelMatrix *_I_star,
                     CompartmentalModelMatrix *_E,
                     InitData *_A0,
                     double *_p_ei,
                     double *_gamma_ei,
                     double _priorAlpha,
                     double _priorBeta,
                     double _sliceWidth)
    {

        context = new ModelContext*;
        I_star = new CompartmentalModelMatrix*;
        E = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_ei = new double*;
        gamma_ei = new double*;
        priorAlpha = new double;
        priorBeta = new double;
        value = new long double;
        samples = new int;
        accepted = new int; 
        sliceWidth = new double;
        varLen = new int;
        *varLen = 1;
        *samples = 0;
        *accepted = 0;
        *gamma_ei = _gamma_ei;
        *context = _context;
        *I_star = _I_star;
        *E = _E;
        *A0 = _A0;
        *p_ei = _p_ei;
        *priorAlpha = _priorAlpha;
        *priorBeta = _priorBeta;
        *value = -1.0;
        *sliceWidth = _sliceWidth;

        // Set up samplers
        samplers = new std::vector<Sampler*>();
        currentSampler = new Sampler*;
        samplers -> push_back(new ParameterSingleMetropolisSampler(*context, this, *gamma_ei));
        samplers -> push_back(new ParameterJointMetropolisSampler(*context, this, *gamma_ei));
        samplers -> push_back(new ParameterJointMetropolisSampler_OCL(*context, this, *gamma_ei));

    }
    FC_Gamma_EI::~FC_Gamma_EI()
    {
        while((samplers -> size()) != 0){delete (*samplers).back(); (*samplers).pop_back();}

        delete samplers;
        delete varLen;
        delete sliceWidth;
        delete I_star;
        delete gamma_ei;
        delete E;
        delete A0;
        delete p_ei;
        delete value;
        delete priorAlpha;
        delete priorBeta;
        delete context;
        delete samples;
        delete accepted;
    }

    double FC_Gamma_EI::evalPrior()
    {
        return((*context) -> random -> dgamma(**gamma_ei, *priorAlpha, 1/(*priorBeta)));
    }

    int FC_Gamma_EI::evalCPU()
    { 
        *value = 0.0;
        int i, j, compIdx;
        int nLoc = *((*E) -> ncol);
        int nTpts = *((*E) -> nrow);
        for (i = 0; i < nLoc; i++)    
        {
            compIdx = i*nTpts;
            for (j = 0; j < nTpts; j++)     
            {
                *value += (*context) -> random -> dbinom(((*I_star) -> data)[compIdx],
                                                         ((*E) -> data)[compIdx],
                                                         (*p_ei)[j]);  
                compIdx++;
            }
        } 
        *value += evalPrior();

        // Catch invalid values, nans etc. 
        if (!std::isfinite(*value))
        {
            *value = -INFINITY;
        }
        return(0);
    }

    int FC_Gamma_EI::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_Gamma_EI::calculateRelevantCompartments()
    {
        (*context) -> calculateP_EI_CPU();
        return(0);
    }
    int FC_Gamma_EI::calculateRelevantCompartments_OCL()
    {
        return(calculateRelevantCompartments());
    }

    void FC_Gamma_EI::sample(int verbose)
    {
        if (verbose){lssCout << "Sampling gamma_EI\n";}
        (*currentSampler) -> drawSample();
    }

    long double FC_Gamma_EI::getValue()
    {
        return(*(this -> value));
    }
    void FC_Gamma_EI::setValue(long double val)
    {
        *(this -> value) = val;
    }
}
