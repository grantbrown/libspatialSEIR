#include<math.h>
#include<cstring>
#include<vector>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FC_P_SE.hpp>
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
    FC_P_SE::FC_P_SE(ModelContext *_context,
                     CompartmentalModelMatrix *_E_star, 
                     CompartmentalModelMatrix *_S, 
                     InitData *_A0,
                     CovariateMatrix *_X,
                     double *_p_se, 
                     double *_beta, 
                     double *_rho,
                     double _sliceWidth,
                     double _priorPrecision,
                     double _priorRhoAlpha,
                     double _priorRhoBeta)
    {
        nBeta = new int;
        nRho = new int;
        varLen = new int;
        *nBeta = (*((_X) -> ncol_x) + *((_X) -> ncol_z));
        *nRho = (_context) -> scaledDistMatrices -> size();
        *varLen = *nBeta + *nRho;

        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        X = new CovariateMatrix*;
        p_se = new double*;
        beta = new double*;
        rho = new double*;
        priorRhoAlpha = new double;
        priorRhoBeta = new double;
        sliceWidth = new double[*varLen];
        priorPrecision = new double;
        value = new long double;
        samples = new int;
        accepted = new int[*varLen]; 
        combinedParams = new double[*varLen];
        *samples = 0;
        memset(accepted, 0, (*nBeta)*sizeof(int)); 
        int i;
        for (i = 0; i < *nBeta; i++)
        {
            sliceWidth[i] = _sliceWidth;       
        }

        *context = _context;
        *E_star = _E_star;
        *S = _S;
        *A0 = _A0;
        *X = _X;
        *p_se = _p_se;
        *beta = _beta;
        *rho = _rho;
        *priorPrecision = _priorPrecision;
        *priorRhoAlpha = _priorRhoAlpha;
        *priorRhoBeta = _priorRhoBeta;
        memcpy(combinedParams, *beta, (*nBeta)*sizeof(double));
        memcpy(&(combinedParams[*nBeta]), *rho, (*nRho)*sizeof(double));
 
        *value = -1.0;

        // Set up samplers
        samplers = new std::vector<Sampler*>();
        currentSampler = new Sampler*;
        samplers -> push_back(new ParameterSingleMetropolisSampler(*context, this, combinedParams));
        samplers -> push_back(new ParameterJointMetropolisSampler(*context, this, combinedParams));
        samplers -> push_back(new ParameterJointMetropolisSampler_OCL(*context, this, combinedParams));
        // Decorrelation sampler is only partial update. 
        samplers -> push_back(new ParameterDecorrelationSampler(*context, this, combinedParams, ((*context) -> X), *nBeta));
    }

    FC_P_SE::~FC_P_SE()
    {
        while((samplers -> size()) != 0){delete (*samplers).back(); (*samplers).pop_back();}
        delete samplers;
        delete combinedParams;
        delete E_star;
        delete S;
        delete A0;
        delete X;
        delete p_se;
        delete beta;
        delete rho;
        delete nBeta;
        delete nRho;
        delete value;
        delete varLen;
        delete[] sliceWidth;
        delete priorPrecision;
        delete priorRhoAlpha;
        delete priorRhoBeta;
        delete context;
        delete samples;
        delete[] accepted;
    }

    double FC_P_SE::evalPrior()
    {
        double out = 0.0;
        double rhoSum = 0.0;
        double rhoVal;
        int i;
        for (i = 0; i < *nBeta; i++)
        {
            out -= pow((*beta)[i],2)*(*priorPrecision); // Generalize to allow different prior precisions. 
        }
        for (i = 0; i < *nRho; i++)
        {
            rhoVal = (*rho)[i];
            out += (rhoVal < 0 ? -INFINITY : ((*context) -> random -> dbeta(rhoVal, *priorRhoAlpha, *priorRhoBeta)));
            rhoSum += rhoVal;
        }
        out += ((rhoSum > 0 && rhoSum < 1 ? 0 : -INFINITY));
        return(out);
    }

    int FC_P_SE::evalCPU()
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

    int FC_P_SE::evalOCL()
    {
        // Not Implemented
        return(evalCPU());
    }
    int FC_P_SE::calculateRelevantCompartments()
    {
        memcpy(combinedParams, *beta, (*nBeta)*sizeof(double));
        memcpy((&(combinedParams[(*nBeta)])), *rho, (*nRho)*sizeof(double));
        ((*context) -> calculateP_SE_CPU());
        return(0);
    }
    int FC_P_SE::calculateRelevantCompartments_OCL()
    {
        memcpy(combinedParams, *beta, (*nBeta)*sizeof(double));
        memcpy((&(combinedParams[(*nBeta)])), *rho, (*nRho)*sizeof(double));
        ((*context) -> calculateP_SE_OCL());
        return(0);

    }

    void FC_P_SE::sample(int verbose)
    {
        if (verbose){lssCout << "Sampling Beta \n";}
        (*currentSampler) -> drawSample();
    }

    long double FC_P_SE::getValue()
    {
        return(*(this -> value));
    }
    void FC_P_SE::setValue(long double val)
    {
        *(this -> value) = val;
    }

}
