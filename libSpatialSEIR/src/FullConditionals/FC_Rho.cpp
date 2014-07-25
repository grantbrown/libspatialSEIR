#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
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
        varLen = new int;
        *varLen = 1;
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

        // Set up samplers
        samplers = new std::vector<Sampler*>();
        currentSampler = new Sampler*;
        samplers -> push_back(new ParameterSingleMetropolisSampler(*context, this, *rho));
        samplers -> push_back(new ParameterJointMetropolisSampler(*context, this, *rho));
        samplers -> push_back(new ParameterJointMetropolisSampler_OCL(*context, this, *rho));

    }
    FC_Rho::~FC_Rho()
    {
        while((samplers -> size()) != 0){delete (*samplers).back(); (*samplers).pop_back();}

        delete samplers;
        delete varLen;
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
        int i, j, Estar_val, S_val, compIdx;
        double p_se_val;
        int nLoc = *((*S) -> ncol);
        int nTpts = *((*S) -> nrow);

        for (i = 0; i < nLoc; i++)    
        {
            compIdx = i*nTpts;
            for (j = 0; j < nTpts; j++)     
            {
                Estar_val = ((*E_star) -> data)[compIdx];
                S_val = ((*S) -> data)[compIdx];
                p_se_val = (*p_se)[compIdx];
                *value += ((*context) -> random -> dbinom(Estar_val, S_val, p_se_val));
                compIdx++;
            }
        } 
        *value += (**rho > 0 && **rho < 1 ? 0 : -INFINITY); // Generalize to allow informative priors. 
                                                        // Prior specification in this area needs work. 
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

    void FC_Rho::sample(int verbose)
    {
        if (verbose){std::cout << "Sampling rho\n";}
        (*currentSampler) -> drawSample();
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
