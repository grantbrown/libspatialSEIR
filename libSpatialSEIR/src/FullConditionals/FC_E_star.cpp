#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FC_E_star.hpp>
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
     * Implement the full conditional distribution for E_Star
     *
     */    
    
    FC_E_Star::FC_E_Star(ModelContext *_context,
                         CompartmentalModelMatrix *_E_star,
                         CompartmentalModelMatrix *_E,  
                         CompartmentalModelMatrix *_S,
                         CompartmentalModelMatrix *_I_star,
                         CovariateMatrix *_X,
                         InitData *_A0,
                         double *_p_se,
                         double *_p_ei,
                         double *_rho,
                         double *_beta,
                         double _steadyStateConstraintPrecision,
                         double _sliceWidth) 
    {

        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        E = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        I_star = new CompartmentalModelMatrix*;
        X = new CovariateMatrix*;
        A0 = new InitData*;
        p_se = new double*;
        p_ei = new double*;
        rho = new double*;
        beta = new double*;
        sliceWidth = new double;
        steadyStateConstraintPrecision = new double;
        value = new long double;
        samples = new int;
        accepted = new int; 
        *samples = 0;
        *accepted = 0;
       
        *context = _context;
        *E_star = _E_star;
        *E = _E;
        *S = _S;
        *I_star = _I_star;
        *X = _X;
        *A0 = _A0;
        *p_se = _p_se;
        *p_ei = _p_ei;
        *rho = _rho;
        *beta = _beta;
        *sliceWidth = _sliceWidth;
        *steadyStateConstraintPrecision = _steadyStateConstraintPrecision;
        *value = -1.0;

        int nvals = *((*context) -> S_star -> nrow) * (*((*context) -> S_star -> ncol));

        samplers = new std::vector<Sampler*>();
        currentSampler = new Sampler*;
        samplers -> push_back(new CompartmentMetropolisSampler(*context, this, (*E_star) -> data));
        samplers -> push_back(new IndexedCompartmentMetropolisSampler(*context, this, (*E_star) -> data));
        samplers -> push_back(new CompartmentMetropolisSampler_OCL(*context, this, (*E_star) -> data));
        samplers -> push_back(new CompartmentBinomialMetropolisSampler(*context, this, (*E_star) -> data, 
                                                                      (*S) -> data, (*E) -> data, *p_se, nvals));
        samplers -> push_back(new CompartmentBinomialSliceSampler(*context, this, (*E_star) -> data, 
                                                                      (*S) -> data, (*E) -> data, *p_se, nvals));

    }

    FC_E_Star::~FC_E_Star()
    {
        while((samplers -> size()) != 0){delete (*samplers).back(); (*samplers).pop_back();}

        delete samplers;
        delete E_star;
        delete E;
        delete S;
        delete I_star;
        delete X;
        delete A0;
        delete p_se;
        delete p_ei;
        delete rho;
        delete beta;
        delete value;
        delete sliceWidth;
        delete steadyStateConstraintPrecision; 
        delete context;
        delete samples;
        delete accepted;
    }

    int FC_E_Star::evalCPU()
    {
        int i,j, compIdx;
        int nTpts = *((*S) -> nrow); 
        int nLoc = *((*S) -> ncol); 
        double p_ei_val; 
        double p_se_val;
        int S_val, E_val, Estar_val, Istar_val;
        long double output = 0.0;
        long unsigned int E_star_sum;
        long unsigned int I_star_sum;
        int64_t aDiff; 

        for (i = 0; i<nLoc; i++)
        {
            compIdx = i*nTpts;
            for (j = 0; j < nTpts; j++)
            {
                Estar_val = ((*E_star) -> data)[compIdx];
                S_val = ((*S) -> data)[compIdx];
                E_val = ((*E) -> data)[compIdx];
                Istar_val = ((*I_star) -> data)[compIdx];
                p_se_val = (*p_se)[compIdx];
                p_ei_val = (*p_ei)[j]; 

                if (Estar_val < 0 || Estar_val > S_val || 
                        Istar_val > E_val)

                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                {
                    output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)) +    
                               ((*context) -> random -> dbinom(Istar_val, E_val, p_ei_val))); 
                }
                compIdx++;
            }

            if (*steadyStateConstraintPrecision > 0)
            {
                E_star_sum = (*E_star)->marginSum(2,i);
                I_star_sum = (*I_star)->marginSum(2,i);
                aDiff = (E_star_sum > I_star_sum ? E_star_sum - I_star_sum : I_star_sum - E_star_sum)/nTpts;
                output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);
            }
        }

        if (!std::isfinite(output))
        {
            *value = -INFINITY;
            return(-1);
        }
        else
        {
            *value = output;
        }
    
       return 0;
    }

    int FC_E_Star::evalCPU(int startLoc, int startTpt)
    {
        int i,j, compIdx;
        int nTpts = *((*S) -> nrow); 
        int nLoc = *((*S) -> ncol); 
        double p_ei_val; 
        double p_se_val;
        int S_val, E_val, Estar_val, Istar_val;
        long double output = 0.0;
        long unsigned int E_star_sum;
        long unsigned int I_star_sum;
        int64_t aDiff; 

        i = startLoc;
        compIdx = i*nTpts + startTpt;
        for (j = startTpt; j < nTpts; j++)
        {
            Estar_val = ((*E_star) -> data)[compIdx];
            S_val = ((*S) -> data)[compIdx];
            E_val = ((*E) -> data)[compIdx];
            Istar_val = ((*I_star) -> data)[compIdx];
            p_se_val = (*p_se)[compIdx];
            p_ei_val = (*p_ei)[j]; 

            if (Estar_val < 0 || Estar_val > S_val || 
                    Istar_val > E_val)

            {
                *value = -INFINITY;
                return(-1);
            }
            else
            {
                output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)) +    
                           ((*context) -> random -> dbinom(Istar_val, E_val, p_ei_val))); 
            }
            compIdx++;
        }

        if (*steadyStateConstraintPrecision > 0)
        {
            E_star_sum = (*E_star)->marginSum(2,i);
            I_star_sum = (*I_star)->marginSum(2,i);
            aDiff = (E_star_sum > I_star_sum ? E_star_sum - I_star_sum : I_star_sum - E_star_sum)/nTpts;
            output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);
        }

        if (!std::isfinite(output))
        {
            *value = -INFINITY;
            return(-1);
        }
        else
        {
            *value = output;
        }
    
       return 0;
    }


    int FC_E_Star::evalOCL()
    {
        //NOT IMPLEMENTED
        return(evalCPU());
    }
    int FC_E_Star::calculateRelevantCompartments()
    {
        (*context) -> calculateE_CPU();
        (*context) -> calculateS_givenE_CPU();
        return(0);
    }
    int FC_E_Star::calculateRelevantCompartments(int i, int j)
    {
        (*context) -> calculateE_CPU(i,j);
        (*context) -> calculateS_givenE_CPU(i,j);
        return(0);
    }
    int FC_E_Star::calculateRelevantCompartments_OCL()
    {
        (*context) -> calculateE_CPU();
        (*context) -> calculateS_givenE_CPU();
        return(0);
    }

    void FC_E_Star::sample(int verbose)
    {
        if (verbose){lssCout << "Sampling E_star\n";}
        (*currentSampler) -> drawSample();
    }

    long double FC_E_Star::getValue()
    {
        return(*(this -> value));
    }
    void FC_E_Star::setValue(long double val)
    {
        *(this -> value) = val;
    }
}
