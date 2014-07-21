#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FC_S0.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;


    /*
     *
     * Implement Full Conditional for S0
     *
     */

    FC_S0::FC_S0(ModelContext* _context,
                 CompartmentalModelMatrix *_S,
                 CompartmentalModelMatrix *_E,
                 CompartmentalModelMatrix *_E_star,
                 CompartmentalModelMatrix *_I_star,
                 InitData *_A0,
                 double *_p_se,
                 double *_p_ei,
                 double _sliceWidth)
    {

        context = new ModelContext*;
        S = new CompartmentalModelMatrix*;
        E = new CompartmentalModelMatrix*;
        E_star = new CompartmentalModelMatrix*;
        I_star = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_se = new double*;
        p_ei = new double*;
        sliceWidth = new double;
        value = new long double;
        samples = new int; *samples = 0;
        accepted = new int; *accepted = 0;

        *context = _context;
        *S = _S;
        *E = _E;
        *E_star = _E_star;
        *I_star = _I_star;
        *A0 = _A0;
        *p_se = _p_se;
        *p_ei = _p_ei;
        *sliceWidth = _sliceWidth;

        // Set up samplers
        samplers = new std::vector<Sampler*>();
        currentSampler = new Sampler*;
        samplers -> push_back(new InitCompartmentMetropolisSampler(*context, this, (*A0) -> S0));
        samplers -> push_back(new IndexedInitCompartmentMetropolisSampler(*context, this, (*A0) -> S0));


    }
    FC_S0::~FC_S0()
    {
        while((samplers -> size()) != 0){delete (*samplers).back(); (*samplers).pop_back();}
        delete samplers;
        delete context;
        delete S;
        delete E;
        delete I_star;
        delete E_star;
        delete p_ei;
        delete p_se;
        delete sliceWidth;
        delete value;
        delete samples;
        delete accepted;
    }
    
    int FC_S0::evalCPU()
    {

        int i,j,compIdx;
        int nLoc = *((*S) -> ncol); 
        int nTpts = *((*S) -> nrow); 
        double p_se_val, p_ei_val;
        int S_val, E_val, Istar_val, Estar_val;
        long double output = 0.0;
        for (i = 0; i<nLoc; i++)
        {
            if (((*A0) -> S0)[i] < 0 || 
                ((*A0) -> E0)[i] < 0)
            {
                *value = -INFINITY;
                return(-1);
            }

            compIdx = i*nTpts;
            for (j = 0; j < nTpts; j++)
            {

                p_ei_val = (*p_ei)[j];       
                S_val = ((*S) -> data)[compIdx];
                E_val = ((*E) -> data)[compIdx];
                Istar_val = ((*I_star)-> data)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                p_se_val = (*p_se)[compIdx];
                if (Estar_val > S_val || Istar_val > E_val)
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



    int FC_S0::evalCPU(int startLoc)
    {

        int j, compIdx;
        int nTpts = *((*S) -> nrow); 
        double p_se_val, p_ei_val;
        int S_val, E_val, Istar_val, Estar_val;
        long double output = 0.0;
        
        if (((*A0) -> S0)[startLoc] < 0 || 
            ((*A0) -> E0)[startLoc] < 0)
        {
            *value = -INFINITY;
            return(-1);
        }

        compIdx = startLoc*nTpts;
        for (j = 0; j < nTpts; j++)
        {
            p_ei_val = (*p_ei)[j];
            S_val = ((*S) -> data)[compIdx];
            E_val = ((*E) -> data)[compIdx];
            Istar_val = ((*I_star)-> data)[compIdx];
            Estar_val = ((*E_star) -> data)[compIdx];
            p_se_val = (*p_se)[compIdx];
            if (Estar_val > S_val || Istar_val > E_val)
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

    int FC_S0::evalOCL()
    {
        // Not Implemented
        evalCPU();
        return(0);
    }

    void FC_S0::sample(int verbose)
    {
        if (verbose){std::cout << "Sampling S0\n";}
        (*currentSampler) -> drawSample();
    }

    long double FC_S0::getValue()
    {
        return(*value);
    }

    void FC_S0::setValue(long double val)
    {
        *(this -> value) = val;
    }

    int FC_S0::calculateRelevantCompartments()
    {
        (*context) -> calculateS_CPU();
        (*context) -> calculateE_givenI_CPU(); 
        return(0);
    }


    int FC_S0::calculateRelevantCompartments_OCL()
    {
        (*context) -> calculateS_CPU();
        (*context) -> calculateE_givenI_CPU(); 
        return(0);
    }

    int FC_S0::calculateRelevantCompartments(int startLoc)
    {
        (*context) -> calculateS_CPU(startLoc,0);
        (*context) -> calculateE_givenI_CPU(startLoc,0);
        return(0);
    }

    void FC_S0::printDebugInfo(int startLoc)
    {
        std::cout << "Error Sampling S0, location: " << startLoc << ", value: " << ((*A0) -> S0)[startLoc] << "\n";
        int j, compIdx;
        int nTpts = *((*S) -> nrow); 
        double p_se_val, p_ei_val;
        int S_val, E_val, Istar_val, Estar_val;
        long double output = 0.0;
        
        if (((*A0) -> S0)[startLoc] < 0 || 
            ((*A0) -> E0)[startLoc] < 0)
        {
            std::cout << "Invalid Value.\n";
            return;
        }

        compIdx = startLoc*nTpts;
        for (j = 0; j < nTpts; j++)
        {
            p_ei_val = (*p_ei)[j];
            S_val = ((*S) -> data)[compIdx];
            E_val = ((*E) -> data)[compIdx];
            Istar_val = ((*I_star)-> data)[compIdx];
            Estar_val = ((*E_star) -> data)[compIdx];
            p_se_val = (*p_se)[compIdx];
            if (Estar_val > S_val || Istar_val > E_val)
            {
                std::cout << "Bounds error detected, time point: " << j << "\n";
                std::cout << "S: " << S_val << "\n";
                std::cout << "E: " << E_val << "\n";
                std::cout << "E_star: " << Estar_val << "\n";
                std::cout << "I_star: " << Istar_val << "\n";
                return;
            }
            else
            {
                output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)) +    
                            ((*context) -> random -> dbinom(Istar_val, E_val, p_ei_val)));

                if (!std::isfinite(output))
                {
                    std::cout << "Calculation Error Detected, time point: " << j << "\n";
                    return;
                }
            }
            compIdx++;
        } 
       return;
    }
}
