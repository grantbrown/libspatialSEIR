#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FC_R0.hpp>
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
     *
     * Implement FC for R0
     * 
     *
     */

    FC_R0::FC_R0(ModelContext* _context,   
                 CompartmentalModelMatrix *_R,
                 CompartmentalModelMatrix *_S,
                 CompartmentalModelMatrix *_S_star,
                 CompartmentalModelMatrix *_E_star,
                 CompartmentalModelMatrix *_R_star,
                 InitData *_A0,
                 double *_p_rs,
                 double *_p_se,
                 double _sliceWidth)
    {
        context = new ModelContext*;
        S = new CompartmentalModelMatrix*;
        R = new CompartmentalModelMatrix*;
        S_star = new CompartmentalModelMatrix*;
        E_star = new CompartmentalModelMatrix*;
        R_star = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_rs = new double*;
        p_se = new double*;
        sliceWidth = new double;
        value = new long double;
        samples = new int; *samples = 0; 
        accepted = new int; *accepted = 0;

        *context = _context;
        *S = _S;
        *R = _R;
        *S_star = _S_star;
        *E_star = _E_star;
        *R_star = _R_star;
        *A0 = _A0;
        *p_se = _p_se;
        *p_rs = _p_rs;
        *sliceWidth = _sliceWidth;

        // Set up samplers
        samplers = new std::vector<Sampler*>();
        currentSampler = new Sampler*;
        samplers -> push_back(new InitCompartmentMetropolisSampler(*context, this, (*A0) -> R0));
        samplers -> push_back(new IndexedInitCompartmentMetropolisSampler(*context, this, (*A0) -> R0));
        samplers -> push_back(new InitCompartmentMetropolisSampler_OCL(*context, this, (*A0) -> R0));
    }
    FC_R0::~FC_R0()
    {
        while((samplers -> size()) != 0){delete (*samplers).back(); (*samplers).pop_back();}
        delete samplers;
        delete context;
        delete S;
        delete R;
        delete S_star;
        delete E_star; 
        delete R_star;
        delete A0;
        delete sliceWidth;
        delete p_se;
        delete p_rs;
        delete value;
        delete samples;
        delete accepted;
    }
    

    int FC_R0::evalCPU()
    {

        int i, j, compIdx;
        int nTpts = *((*R) -> nrow);
        int nLoc = *((*R) -> ncol);


        long double output = 0.0;
        
        double p_se_val;
        double p_rs_val;
        int Sstar_val, Estar_val, R_val, S_val;   

        // Is p_rs meaningful?
        if ((*context) -> config -> reinfectionMode <= 2)
        {
            for (i = 0; i < nLoc;i++)
            {
                if (((*A0) -> S0)[i] < 0 || 
                    ((*A0) -> R0)[i] < 0)
                {
                    *value = -INFINITY;
                    return(-1);
                }

                compIdx = i*nTpts;

                for (j = 0; j < nTpts; j++)
                {
                    Estar_val = ((*E_star) -> data)[compIdx];
                    Sstar_val = ((*S_star)->data)[compIdx];
                    R_val = ((*R) ->data)[compIdx];
                    S_val = ((*S) ->data)[compIdx];
                    p_rs_val = (*p_rs)[j];

                    if (Estar_val > S_val || 
                            Sstar_val > R_val)
                    {
                        *value = -INFINITY;
                        return(-1);
                    }
                    else
                    {
                        output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)) + 
                                   ((*context) -> random -> dbinom(Sstar_val, R_val, p_rs_val)));

                    }
                    compIdx++;
                } 
            }
        }
        else 
        {
            for (i = 0; i < nLoc;i++)
            {

                if (((*A0) -> S0)[i] < 0 || 
                    ((*A0) -> R0)[i] < 0)
                {
                    *value = -INFINITY;
                    return(-1);
                }

                compIdx = i*nTpts;

                for (j = 0; j < nTpts; j++)
                {
                    Estar_val = ((*E_star) -> data)[compIdx];
                    R_val = ((*R) ->data)[compIdx];
                    S_val = ((*S) ->data)[compIdx];
                    p_rs_val = (*p_rs)[j];

                    if (Estar_val > S_val)
                    {
                        *value = -INFINITY;
                        return(-1);
                    }
                    else
                    {
                        output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)));

                    }
                    compIdx++;
                } 
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

    int FC_R0::evalCPU(int startLoc)
    {

        int j, compIdx;
        int nTpts = *((*R) -> nrow);

        long double output = 0.0;
        
        double p_se_val;
        double p_rs_val;
        int Sstar_val, Estar_val, R_val, S_val;   

        if (((*A0) -> S0)[startLoc] < 0 || 
            ((*A0) -> R0)[startLoc] < 0)
        {
            *value = -INFINITY;
            return(-1);
        }

        compIdx = startLoc*nTpts;

        // Is p_rs meaningful?
        if ((*context) -> config -> reinfectionMode <= 2)
        {
            for (j = 0; j < nTpts; j++)
            {
                Estar_val = ((*E_star) -> data)[compIdx];
                Sstar_val = ((*S_star)->data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                S_val = ((*S) ->data)[compIdx];
                p_rs_val = (*p_rs)[j];

                if (Estar_val > S_val || 
                        Sstar_val > R_val)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                {
                    output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)) + 
                               ((*context) -> random -> dbinom(Sstar_val, R_val, p_rs_val)));

                }
                compIdx++;
            } 
        }
        else 
        {
            for (j = 0; j < nTpts; j++)
            {
                Estar_val = ((*E_star) -> data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                S_val = ((*S) ->data)[compIdx];
                p_rs_val = (*p_rs)[j];

                if (Estar_val > S_val)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                {
                    output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)));

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

    int FC_R0::evalOCL()
    {
        // Not Implemented
        return(evalCPU());
    }

    void FC_R0::sample(int verbose)
    {
        if (verbose){lssCout << "Sampling R0\n";}
        (*currentSampler) -> drawSample();
    }

    long double FC_R0::getValue()
    {
        return(*value);
    }

    void FC_R0::setValue(long double val)
    {
        *(this -> value) = val;
    }

    int FC_R0::calculateRelevantCompartments()
    {
        (*context) -> calculateR_CPU();
        (*context) -> calculateS_givenE_CPU();
        return(0);
    }
    int FC_R0::calculateRelevantCompartments_OCL()
    {
        (*context) -> calculateR_CPU();
        (*context) -> calculateS_givenE_CPU();
        return(0);
    }
    int FC_R0::calculateRelevantCompartments(int startLoc)
    {
        (*context) -> calculateR_CPU(startLoc, 0);
        (*context) -> calculateS_givenE_CPU(startLoc,0);
        return(0);
    }

    void FC_R0::printDebugInfo(int loc)
    {
        lssCout << "Error Sampling R0, location: " << loc << ", value: " << ((*A0) -> R0)[loc] << "\n";
    }
}
