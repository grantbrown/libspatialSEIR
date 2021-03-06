#include<math.h>
#include<cstring>
#include<vector>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FC_Beta.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>
#include<IOProvider.hpp>


namespace SpatialSEIR
{
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
                     double *_rho,
                     double _sliceWidth,
                     double *_priorMean,
                     double *_priorPrecision)
    {
        int nBeta = (*((_X) -> ncol_x));
        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        X = new CovariateMatrix*;
        p_se = new double*;
        beta = new double*;
        rho = new double*;
        sliceWidth = new double[nBeta];
        priorPrecision = new double[nBeta];
        priorMean = new double[nBeta];
        value = new long double;
        samples = new int;
        accepted = new int[nBeta]; 
        varLen = new int;
        *varLen = nBeta;
        *samples = 0;
        memset(accepted, 0, nBeta*sizeof(int)); 
        int i;
        for (i = 0; i < nBeta; i++)
        {
            sliceWidth[i] = _sliceWidth;       
            priorPrecision[i] = _priorPrecision[i];
            priorMean[i] = _priorMean[i];
        }

        *context = _context;
        *E_star = _E_star;
        *S = _S;
        *A0 = _A0;
        *X = _X;
        *p_se = _p_se;
        *beta = _beta;
        *rho = _rho;
        *value = -1.0;

        // Set up samplers
        samplers = new std::vector<Sampler*>();
        currentSampler = new Sampler*;
        samplers -> push_back(new ParameterSingleMetropolisSampler(*context, this, *beta));
        samplers -> push_back(new ParameterJointMetropolisSampler(*context, this, *beta));
        samplers -> push_back(new ParameterJointMetropolisSampler_OCL(*context, this, *beta));
        samplers -> push_back(new ParameterDecorrelationSampler(*context, this, *beta, (*context) -> X));
        samplers -> push_back(new ParameterNullSampler());
    }

    FC_Beta::~FC_Beta()
    {
        while((samplers -> size()) != 0){delete (*samplers).back(); (*samplers).pop_back();}
        delete samplers;
        delete E_star;
        delete S;
        delete A0;
        delete X;
        delete p_se;
        delete beta;
        delete rho;
        delete value;
        delete varLen;
        delete[] sliceWidth;
        delete[] priorPrecision;
        delete[] priorMean;
        delete context;
        delete samples;
        delete[] accepted;
    }

    double FC_Beta::evalPrior()
    {
        double out = 0.0;
        int i;
        for (i = 0; i < (*((*X) -> ncol_x)); i++)
        {
            out -= pow((*beta)[i] - priorMean[i],2)*(priorPrecision[i]);
        }
        return(out);
    }

    int FC_Beta::evalCPU()
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*S) -> ncol);
        int nTpts = *((*S) -> nrow);
        double term1, term2;
        term1 = 0.0; term2 = 0.0;

        for (i = 0; i < nLoc; i++)    
        {
            compIdx = i*nTpts;
            for (j = 0; j < nTpts; j++)     
            {
                // todo: clean up
                tmp = ((*E_star) -> data)[compIdx];
                term1 += (*context) -> random -> dbinom(tmp, ((*S) -> data)[compIdx], (*p_se)[compIdx]);
                //term1 += std::log((*p_se)[compIdx])*tmp; 
                //term2 += std::log(1-(*p_se)[compIdx])*(((*S) -> data)[compIdx] - tmp);
                compIdx++;
            }
        } 

        *value = term1 + term2 + evalPrior();
        // Catch invalid values, nans etc. 
        if (!std::isfinite(*value))
        {
            *value = -INFINITY;
        }
        return(0);
    }

    int FC_Beta::evalOCL()
    {
        // Not Implemented
        return(evalCPU());
    }
    int FC_Beta::calculateRelevantCompartments()
    {
        ((*context) -> calculateP_SE_CPU());
        return(0);
    }
    int FC_Beta::calculateRelevantCompartments_OCL()
    {
        ((*context) -> calculateP_SE_OCL());
        return(0);

    }

    void FC_Beta::sample(int verbose)
    {
        if (verbose){lssCout << "Sampling Beta \n";}
        (*currentSampler) -> drawSample();
    }

    long double FC_Beta::getValue()
    {
        return(*(this -> value));
    }
    void FC_Beta::setValue(long double val)
    {
        *(this -> value) = val;
    }

}
