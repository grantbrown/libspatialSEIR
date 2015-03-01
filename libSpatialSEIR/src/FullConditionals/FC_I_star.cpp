#include<math.h>
#include<cstring>
#include<vector>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FC_I_star.hpp>
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
     * Implement the full conditional distribution for I_Star
     *
     */    
    FC_I_Star::FC_I_Star(ModelContext *_context,
                         CompartmentalModelMatrix *_E_star,
                         CompartmentalModelMatrix *_I_star,
                         CompartmentalModelMatrix *_R_star,
                         CompartmentalModelMatrix *_E,
                         CompartmentalModelMatrix *_I,                       
                         CompartmentalModelMatrix *_S,                       
                         InitData *_A0,
                         double *_p_ei,
                         double *_p_ir,
                         double *_p_se,
                         double _steadyStateConstraintPrecision,
                         double _sliceWidth)
    {

        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        I_star = new CompartmentalModelMatrix*;
        R_star = new CompartmentalModelMatrix*;
        E = new CompartmentalModelMatrix*;
        I = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;

        A0 = new InitData*;
        p_ei = new double*;
        p_ir = new double*;
        p_se = new double*;
        sliceWidth = new double;
        steadyStateConstraintPrecision = new double;
        value = new long double;
        samples = new int;
        accepted = new int; 
        *samples = 0;
        *accepted = 0;
        *context = _context;
        *E_star = _E_star;
        *I_star = _I_star;
        *R_star = _R_star;
        *E = _E;
        *I = _I;
        *S = _S;

        *A0 = _A0;
        *p_ei = _p_ei;
        *p_ir = _p_ir;
        *p_se = _p_se;
        *sliceWidth = _sliceWidth;
        *steadyStateConstraintPrecision = _steadyStateConstraintPrecision;
        *value = -1.0;

        // Set up samplers
        samplers = new std::vector<Sampler*>();
        currentSampler = new Sampler*;
        samplers -> push_back(new CompartmentMetropolisSampler(*context, this, (*I_star) -> data));
        samplers -> push_back(new IndexedCompartmentMetropolisSampler(*context, this, (*I_star) -> data));
        samplers -> push_back(new CompartmentMetropolisSampler_OCL(*context, this, (*I_star) -> data));
        samplers -> push_back(new CompartmentBinomialMetropolisSampler(*context, this, (*I_star) -> data, (*E) -> data, (*I) -> data, *p_ei, 
                    *((*context) -> S_star -> nrow)));
        samplers -> push_back(new IndexedCompartmentBinomialMetropolisSampler(*context, this, (*I_star) -> data, (*E) -> data, (*I) -> data, *p_ei, 
                    *((*context) -> S_star -> nrow)));
        samplers -> push_back(new CompartmentBinomialMixedSampler(*context, this, (*I_star) -> data, (*E) -> data, (*I) -> data, *p_ei, 
                    *((*context) -> S_star -> nrow)));
        samplers -> push_back(new CompartmentBinomialSliceSampler(*context, this, (*I_star) -> data, (*E) -> data, (*I) -> data, *p_ei, 
                    *((*context) -> S_star -> nrow)));


    }
    FC_I_Star::~FC_I_Star()
    {
        while((samplers -> size()) != 0){delete (*samplers).back(); (*samplers).pop_back();}

        delete samplers;
        delete R_star;
        delete E;
        delete I;
        delete S;
        delete E_star;
        delete I_star;
        delete A0;
        delete p_ei;
        delete p_se;
        delete p_ir;
        delete value;
        delete sliceWidth;
        delete steadyStateConstraintPrecision;
        delete context;
        delete samples;
        delete accepted;
    }

    int FC_I_Star::evalCPU()
    {
        int i,j, compIdx;
        int nTpts = *((*I) -> nrow);
        int nLoc = *((*I) -> ncol);

        long double output = 0.0;
        
        double p_ei_val;
        double p_ir_val;
        double p_se_val;
        int Rstar_val, Istar_val, E_val, I_val, S_val, Estar_val;   
        long unsigned int I_star_sum;
        long unsigned int E_star_sum;
        int64_t aDiff; 

        for (i = 0; i < nLoc; i++)
        {
            compIdx = i*nTpts;
            for (j = 0; j < nTpts; j++)
            {
                Rstar_val = ((*R_star) -> data)[compIdx];
                Istar_val = ((*I_star)->data)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                S_val = ((*S) ->data)[compIdx];
                E_val = ((*E) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];
                p_se_val = (*p_se)[compIdx];
                p_ei_val = (*p_ei)[j];
                p_ir_val = (*p_ir)[j];

                if (Rstar_val < 0 || Istar_val > E_val || 
                        Rstar_val > I_val)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                {
                    output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)) + 
                               ((*context) -> random -> dbinom(Istar_val, E_val, p_ei_val)) + 
                               ((*context) -> random -> dbinom(Rstar_val, I_val, p_ir_val)));
                }
                compIdx++;
            } 
        }

        if (*steadyStateConstraintPrecision > 0)
        {
            I_star_sum = (*I_star)->marginSum(3,-1);
            E_star_sum = (*E_star)->marginSum(3,-1);
            aDiff = (I_star_sum > E_star_sum ? I_star_sum - E_star_sum : E_star_sum - I_star_sum)/(nTpts*nLoc);
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

    int FC_I_Star::evalCPU(int startLoc, int startTpt)
    {
        int i,j, compIdx;
        int nTpts = *((*I) -> nrow);
        int nLoc = *((*I) -> ncol);

        long double output = 0.0;
        
        double p_ei_val;
        double p_ir_val;
        double p_se_val;
        int Istar_val, Estar_val, Rstar_val, E_val, I_val, S_val;   
        long unsigned int I_star_sum;
        long unsigned int E_star_sum;
        int64_t aDiff; 

        i = startLoc;
        compIdx = i*nTpts + startTpt;
        for (j = startTpt; j < nTpts; j++)
        {
            Rstar_val = ((*R_star) -> data)[compIdx];
            Istar_val = ((*I_star)->data)[compIdx];
            E_val = ((*E) ->data)[compIdx];
            I_val = ((*I) ->data)[compIdx];
            p_ei_val = (*p_ei)[j];
            p_ir_val = (*p_ir)[j];

            if (Rstar_val < 0 || Istar_val > E_val || 
                    Rstar_val > I_val)
            {
                *value = -INFINITY;
                return(-1);
            }
            else
            {
                output += (((*context) -> random -> dbinom(Rstar_val, I_val, p_ir_val)) + 
                           ((*context) -> random -> dbinom(Istar_val, E_val, p_ei_val)));
            }
            compIdx++;
        } 

        for (i = 0; i < nLoc; i++)
        {
            compIdx = i*nTpts + startTpt;
            for (j = startTpt; j < nTpts; j++)
            {

                Estar_val = ((*E_star)->data)[compIdx];
                S_val = ((*S) ->data)[compIdx];
                p_se_val = (*p_se)[compIdx];
                output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)));
                compIdx++;
            } 
        }
        if (*steadyStateConstraintPrecision > 0)
        {
            I_star_sum = (*I_star)->marginSum(3,-1);
            E_star_sum = (*E_star)->marginSum(3,-1);
            aDiff = (I_star_sum > E_star_sum ? I_star_sum - E_star_sum : E_star_sum - I_star_sum)/(nTpts*nLoc);
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



    int FC_I_Star::evalOCL()
    {
        // Not yet implemented
        return(evalCPU()); 
    }
    int FC_I_Star::calculateRelevantCompartments()
    {        
        (*context) -> calculateI_CPU();
        (*context) -> calculateE_givenI_CPU();
        ((*context) -> calculateP_SE_CPU());
        return(0);
    }
    int FC_I_Star::calculateRelevantCompartments(int i, int j)
    {        
        (*context) -> calculateI_CPU(i,j);
        (*context) -> calculateE_givenI_CPU(i,j);
        ((*context) -> calculateP_SE_CPU(i,j));
        return(0);
    }


    int FC_I_Star::calculateRelevantCompartments_OCL()
    {
        (*context) -> calculateI_CPU();
        (*context) -> calculateE_givenI_CPU();
        ((*context) -> calculateP_SE_OCL());
        return(0);
    }

    void FC_I_Star::sample(int verbose)
    {
        if (verbose){lssCout << "Sampling I_star\n";}
        (*context) -> cacheP_SE_Calculation();
        (*currentSampler) -> drawSample();
    }

    long double FC_I_Star::getValue()
    {
        return(*(this -> value));
    }
    void FC_I_Star::setValue(long double val)
    {
        *(this -> value) = val;
    }

}
