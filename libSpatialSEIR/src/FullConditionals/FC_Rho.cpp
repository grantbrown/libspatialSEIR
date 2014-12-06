#include<math.h>
#include<cstring>
#include<vector>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FC_Rho.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>
#include<IOProvider.hpp>

namespace SpatialSEIR
{
    FC_Rho::FC_Rho(ModelContext *_context,
                   CompartmentalModelMatrix *_E_star,  
                   CompartmentalModelMatrix *_S,
                   InitData *_A0,
                   CovariateMatrix *_X,
                   double *_p_se,
                   double *_beta,
                   double *_rho,
                   double _sliceWidth,
                   double priorAlpha_rho_,
                   double priorBeta_rho_)
    {
        context = new ModelContext*;
        *context = _context;
        varLen = new int;
        *varLen = (*context) -> scaledDistMatrices -> size(); 
        E_star = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        X = new CovariateMatrix*;
        p_se = new double*;
        beta = new double*;
        rho = new double*;
        sliceWidth = new double[*varLen];
        value = new long double;
        samples = new int;
        accepted = new int[*varLen]; 
        priorAlpha = new double;
        priorBeta = new double;
        *priorAlpha = priorAlpha_rho_;
        *priorBeta = priorBeta_rho_;
        *samples = 0;
     

        *E_star = _E_star;
        *S = _S;
        *A0 = _A0;
        *X = _X;
        *p_se = _p_se;
        *beta = _beta;
        *rho = _rho;
        int i;
        for (i = 0; i < *varLen; i++)
        {
            sliceWidth[i] = _sliceWidth;
            accepted[i] = 0;
        }
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
        delete priorAlpha;
        delete priorBeta;
    }

    double FC_Rho::evalPrior()
    {
        double rhoSum = 0.0;
        double out = 0.0;
        double rhoVal;
        int i;
        for (i = 0; i < *varLen; i++)
        {
            rhoVal = (*rho)[i];
            if (rhoVal < 0){return(-INFINITY);}
            rhoSum += rhoVal; 
            out += (*context) -> random -> dbeta(rhoVal, *priorAlpha, *priorBeta);
        }
         return(rhoSum > 0 && rhoSum < 1 ? out : -INFINITY); 
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
        *value += evalPrior();
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
        if (verbose){lssCout << "Sampling rho\n";}
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
