#include<math.h>
#include<cstring>
#include<vector>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FC_Phi.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>
#include<IOProvider.hpp>

namespace SpatialSEIR
{
    FC_Phi::FC_Phi(ModelContext *_context,
                   /*CompartmentalModelMatrix *_I_star, */ 
                   CompartmentalModelMatrix *_Compartment,  
                   double* _phi,
                   double _priorAlpha,
                   double _priorBeta,
                   int* _Y,
                   double _sliceWidth)
    {
        context = new ModelContext*;
        Compartment = new CompartmentalModelMatrix*;
        phi = new double*;
        priorAlpha = new double;
        priorBeta = new double;
        sliceWidth = new double;
        Y = new int*;
        value = new long double;
        samples = new int;
        accepted = new int; 
        varLen = new int;
        *varLen = 1;
        *samples = 0;
        *accepted = 0;


        *context = _context;
        *Compartment = _Compartment;
        *phi = _phi;
        *priorAlpha = _priorAlpha;
        *priorBeta = _priorBeta;
        *Y = _Y;
        *sliceWidth = _sliceWidth;    
        *value = -1.0;

        // Set up samplers
        samplers = new std::vector<Sampler*>();
        currentSampler = new Sampler*;
        samplers -> push_back(new ParameterSingleMetropolisSampler(*context, this, *phi));
        samplers -> push_back(new ParameterNullSampler());
        samplers -> push_back(new ParameterJointMetropolisSampler(*context, this, *phi));
        samplers -> push_back(new ParameterJointMetropolisSampler_OCL(*context, this, *phi));

    }
    FC_Phi::~FC_Phi()
    {
        while((samplers -> size()) != 0){delete (*samplers).back(); (*samplers).pop_back();}
        delete samplers;
        delete varLen;
        delete Compartment;
        delete phi;
        delete priorAlpha;
        delete priorBeta;
        delete sliceWidth;
        delete value;
        delete context;
        delete samples;
        delete accepted;
    }

    double* FC_Phi::minimumValue()
    {
        // Not Implemented
        return(new double);
    }
    double* FC_Phi::maximumValue()
    {
        // Not Implemented
        return(new double);
    }


    double FC_Phi::evalPrior()
    {
        return((*context) -> random -> dgamma(**phi, *priorAlpha, 1/(*priorBeta)));
    }

    int FC_Phi::evalCPU()
    {
        *value = 0.0;
        int nLoc = *((*Compartment) -> ncol);
        int nTpts = *((*Compartment) -> nrow);
        int maxIdx = nLoc*nTpts;
        int i;
        double phi_val = **phi;
        for (i = 0; i < maxIdx; i++)    
        {
            *value -= 0.5*std::pow((((*Compartment)->data)[i] - (*Y)[i])*(phi_val), 2);
        }
        *value += evalPrior();
        // Catch invalid values, nans etc. 
        if (!std::isfinite(*value))
        {
            *value = -INFINITY;
        }

        return(0);
    }

    int FC_Phi::evalOCL()
    {
        //NOT IMPLEMENTED
        return(evalCPU());
    }
    int FC_Phi::calculateRelevantCompartments()
    {
       return(0); 
    }
    int FC_Phi::calculateRelevantCompartments_OCL()
    {
       return(0); 
    }

    void FC_Phi::sample(int verbose)
    {
        if (verbose){lssCout << "Sampling phi\n";}
        (*currentSampler) -> drawSample();
    }

    long double FC_Phi::getValue()
    {
        return(*(this -> value));
    }
    void FC_Phi::setValue(long double val)
    {
        *(this -> value) = val;
    }
}
