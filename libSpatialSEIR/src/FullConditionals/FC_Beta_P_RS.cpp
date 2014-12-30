#include<math.h>
#include<cstring>
#include<vector>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FC_Beta_P_RS.hpp>
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

        int nBeta = (*((_X) -> ncol_x));

        context = new ModelContext*;
        S_star = new CompartmentalModelMatrix*;
        R = new CompartmentalModelMatrix*;
        X = new CovariateMatrix*;
        A0 = new InitData*;
        p_rs = new double*;
        beta_p_rs = new double*;
        tausq = new double;
        sliceWidth = new double[nBeta];
        value = new long double;
        samples = new int; 
        accepted = new int[nBeta];
        varLen = new int;

        *varLen = nBeta;
        int i;
        for (i = 0; i < nBeta; i++)
        {
            sliceWidth[i] = _sliceWidth;       
        }
        *samples = 0;
        memset(accepted, 0, nBeta*sizeof(int)); 

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

        // Set up samplers
        samplers = new std::vector<Sampler*>();
        currentSampler = new Sampler*;
        samplers -> push_back(new ParameterSingleMetropolisSampler(*context, this, *beta_p_rs));
        samplers -> push_back(new ParameterJointMetropolisSampler(*context, this, *beta_p_rs));
        samplers -> push_back(new ParameterJointMetropolisSampler_OCL(*context, this, *beta_p_rs));
        samplers -> push_back(new ParameterDecorrelationSampler(*context, this, *beta_p_rs, (*context) -> X_pRS));
        samplers -> push_back(new ParameterNullSampler());
    }
    FC_Beta_P_RS::~FC_Beta_P_RS()
    {
        while((samplers -> size()) != 0){delete (*samplers).back(); (*samplers).pop_back();}

        delete samplers;
        delete S_star;
        delete varLen;
        delete R;
        delete X;
        delete beta_p_rs;
        delete tausq;
        delete A0;
        delete p_rs;
        delete value;
        delete[] sliceWidth;
        delete context;
        delete samples;
        delete[] accepted;
    }

    double FC_Beta_P_RS::evalPrior()
    {
        double out = 0.0;
        int nbeta = *((*X) -> ncol_x);
        int j;
        for (j = 0; j < nbeta; j++)
        {
            out -= ((*tausq)/2)*pow((*beta_p_rs)[j],2);
        }
        return(out);
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

    void FC_Beta_P_RS::sample(int verbose)
    {
        if (verbose){lssCout << "Sampling Beta_P_RS\n";}
        (*currentSampler) -> drawSample();
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
