#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FC_I0.hpp>
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
     * Implement full conditional for I0
     *
     */
    
    FC_I0::FC_I0(ModelContext* _context, 
                 CompartmentalModelMatrix *_S,
                 CompartmentalModelMatrix *_I,
                 CompartmentalModelMatrix *_R,
                 CompartmentalModelMatrix *_S_star,
                 CompartmentalModelMatrix *_E_star,
                 CompartmentalModelMatrix *_R_star,
                 InitData *_A0,
                 double *_p_ir,
                 double *_p_rs,
                 double *_p_se,
                 double _sliceWidth)
    {
        context = new ModelContext*;
        S = new CompartmentalModelMatrix*;
        I = new CompartmentalModelMatrix*;
        R = new CompartmentalModelMatrix*;
        S_star = new CompartmentalModelMatrix*;
        E_star = new CompartmentalModelMatrix*;
        R_star = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_ir = new double*;
        p_rs = new double*;
        p_se = new double*;
        sliceWidth = new double;
        value = new long double;
        samples = new int; *samples =0;
        accepted = new int; *accepted = 0;

        *context = _context;
        *S = _S;
        *I = _I;
        *R = _R;
        *S_star = _S_star;
        *E_star = _E_star;
        *R_star = _R_star;
        *A0 = _A0;
        *p_ir = _p_ir;
        *p_se = _p_se;
        *p_rs = _p_rs;
        *sliceWidth = _sliceWidth;

        // Set up samplers
        samplers = new std::vector<Sampler*>();
        currentSampler = new Sampler*;
        samplers -> push_back(new InitCompartmentMetropolisSampler(*context, this, (*A0) -> I0));
        samplers -> push_back(new IndexedInitCompartmentMetropolisSampler(*context, this, (*A0) -> I0));
        samplers -> push_back(new InitCompartmentMetropolisSampler_OCL(*context, this, (*A0)-> I0));

    }
    FC_I0::~FC_I0()
    {
        while((samplers -> size()) != 0){delete (*samplers).back(); (*samplers).pop_back();}

        delete samplers;
        delete context;
        delete S;
        delete I;
        delete R;
        delete S_star;
        delete E_star; 
        delete R_star;
        delete A0;
        delete sliceWidth;
        delete p_ir;
        delete p_se;
        delete p_rs;
        delete value;
        delete samples;
        delete accepted;
    }
    
    int FC_I0::evalCPU()
    {

        int i,j, compIdx;
        int nTpts = *((*R) -> nrow);
        int nLoc = *((*R) -> ncol);

        long double output = 0.0;
        
        double p_se_val;
        double p_rs_val;
        double ln_1m_p_ir;
        int Rstar_val, Sstar_val, Estar_val, R_val, I_val, S_val;   

        // Cache ln_1m_p_ir to prevent work duplication
        // Need to revert this before exit.
        for (j = 0; j<nTpts; j++)
        {
            (*p_ir)[j] = std::log(1-(*p_ir)[j]);
        }

        // Is p_rs meaningful?
        if ((*context) -> config -> reinfectionMode <= 2)
        {
            for (i = 0; i<nLoc; i++)
            {
                if (((*A0) -> I0)[i] < 0 || 
                    ((*A0) -> R0)[i] < 0)
                {
                    *value = -INFINITY;
                    (*context) -> calculateP_IR_CPU();
                    return(-1);
                }
                compIdx = i*nTpts;
                for (j = 0; j < nTpts; j++)
                {
                    ln_1m_p_ir = (*p_ir)[j]; 
                    Estar_val = ((*E_star) -> data)[compIdx];
                    Rstar_val = ((*R_star) -> data)[compIdx];
                    Sstar_val = ((*S_star)->data)[compIdx];
                    S_val = ((*S)->data)[compIdx];
                    R_val = ((*R) ->data)[compIdx];
                    I_val = ((*I) ->data)[compIdx];
                    p_rs_val = (*p_rs)[j];
                    p_se_val = (*p_se)[compIdx];


                    if (Rstar_val > I_val || 
                            Sstar_val > R_val || 
                            p_se_val > 1 || p_se_val < 0)
                    {
                        *value = -INFINITY;
                        (*context) -> calculateP_IR_CPU();
                        return(-1);
                    }
                    else
                    {
                        output +=  (ln_1m_p_ir*(I_val) +
                                    ((*context) -> random -> dbinom(Sstar_val, R_val, p_rs_val)) + 
                                    ((*context) -> random -> choose(I_val, Rstar_val)) + 
                                    ((*context) -> random -> dbinom(Estar_val,S_val, p_se_val)));
                    }
                    compIdx++;
                } 
            }
        }
        else
        {
            for (i = 0; i<nLoc; i++)
            {
                compIdx = i*nTpts;
                if (((*A0) -> I0)[i] < 0 || 
                    ((*A0) -> R0)[i] < 0)
                {
                    *value = -INFINITY;
                    (*context) -> calculateP_IR_CPU();
                    return(-1);
                }

                for (j = 0; j < nTpts; j++)
                {
                    Rstar_val = ((*R_star) -> data)[compIdx];
                    S_val = ((*S)->data)[compIdx];
                    I_val = ((*I) ->data)[compIdx];
                    R_val = ((*R) ->data)[compIdx];
                    p_se_val = (*p_se)[compIdx];
                    Estar_val = ((*E_star) -> data)[compIdx];
                    ln_1m_p_ir = (*p_ir)[j]; 


                    if (Rstar_val > I_val || p_se_val > 1 || p_se_val < 0)
                    {
                        *value = -INFINITY;
                        (*context) -> calculateP_IR_CPU();
                        return(-1);
                    }
                    else
                    {
                        output +=  (ln_1m_p_ir*(I_val) +
                                    ((*context) -> random -> choose(I_val, Rstar_val)) + 
                                    ((*context) -> random -> dbinom(Estar_val,S_val, p_se_val)));
                    }
                    compIdx++;
                } 
            }
        }
        (*context) -> calculateP_IR_CPU();
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


    int FC_I0::evalCPU(int startLoc)
    {
        int i,j, compIdx;
        int nTpts = *((*R) -> nrow);
        int nLoc = *((*R) -> ncol);

        long double output = 0.0;
        
        double p_se_val;
        double p_rs_val;
        double ln_1m_p_ir;
        int Rstar_val, Sstar_val, Estar_val, R_val, I_val, S_val;   
 
        if (((*A0) -> I0)[startLoc] < 0 || 
            ((*A0) -> R0)[startLoc] < 0)
        {
            *value = -INFINITY;
            return(-1);
        }

        // Cache ln_1m_p_ir to prevent work duplication
        // Need to revert this before exit.
        for (j = 0; j<nTpts; j++)
        {
            (*p_ir)[j] = std::log(1-(*p_ir)[j]);
        }

        compIdx = startLoc*nTpts;

        // Is p_rs meaningful?
        if ((*context) -> config -> reinfectionMode <= 2)
        {
            for (j = 0; j < nTpts; j++)
            {
                Rstar_val = ((*R_star) -> data)[compIdx];
                Sstar_val = ((*S_star)->data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];
                p_rs_val = (*p_rs)[j];
                ln_1m_p_ir = (*p_ir)[j];
                if (Rstar_val > I_val || 
                        Sstar_val > R_val)
                {
                    *value = -INFINITY;
                    (*context) -> calculateP_IR_CPU();
                    return(-1);
                }
                else
                {
                    output +=  (ln_1m_p_ir*(I_val) +
                                ((*context) -> random -> dbinom(Sstar_val, R_val, p_rs_val)) + 
                                ((*context) -> random -> choose(I_val, Rstar_val)));
                }
                compIdx++;
            } 
        }
        else
        {
            for (j = 0; j < nTpts; j++)
            {
                Rstar_val = ((*R_star) -> data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];
                ln_1m_p_ir = (*p_ir)[j];

                if (Rstar_val > I_val)
                {
                    *value = -INFINITY;
                    (*context) -> calculateP_IR_CPU();
                    return(-1);
                }
                else
                {
                    output +=  (ln_1m_p_ir*(I_val) +
                                ((*context) -> random -> choose(I_val, Rstar_val)));
                }
                compIdx++;
            } 
        }

        // p_se changes, so need to look at p_se component for all locations and 
        // time points after 0
        for (i = 0; i < nLoc; i++)
        {
            compIdx = i*nTpts;
            for (j = 0; j< nTpts; j++)
            {

                p_se_val = (*p_se)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                if (p_se_val > 1 || p_se_val < 0)
                {
                    *value = -INFINITY;
                    (*context) -> calculateP_IR_CPU();
                    return(-1);
                }

                output += (*context) -> random -> dbinom(Estar_val,S_val, p_se_val);
                compIdx ++; 
            }
        }
        (*context) -> calculateP_IR_CPU();
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

    int FC_I0::evalOCL()
    {
        // Not Implemented
        return(evalCPU());
    }

    void FC_I0::sample(int verbose)
    {
        if (verbose){lssCout << "Sampling I0\n";}
        (*currentSampler) -> drawSample();
    }

    long double FC_I0::getValue()
    {
        return(*value);
    }

    void FC_I0::setValue(long double val)
    {
        *(this -> value) = val;
    }

    int FC_I0::calculateRelevantCompartments()
    {
        (*context) -> calculateI_CPU();
        (*context) -> calculateR_givenS_CPU();
        (*context) -> calculateP_SE_CPU();
        return(0);
    }

    int FC_I0::calculateRelevantCompartments_OCL()
    {
        (*context) -> calculateI_CPU();
        (*context) -> calculateR_givenS_CPU();
        (*context) -> calculateP_SE_OCL();
        return(0);
    }


    int FC_I0::calculateRelevantCompartments(int startLoc)
    {
        (*context) -> calculateI_CPU(startLoc, 0);
        (*context) -> calculateR_givenS_CPU(startLoc, 0);
        (*context) -> calculateP_SE_CPU(startLoc, 0);
        return(0);
    }

    void FC_I0::printDebugInfo(int startLoc)
    {
        lssCout << "Error Sampling I0, location: " << startLoc << ", value: " << ((*A0) -> I0)[startLoc] << "\n";

        int i,j, compIdx;
        int nTpts = *((*R) -> nrow);
        int nLoc = *((*R) -> ncol);

        long double output = 0.0;
        
        double p_se_val;
        double p_rs_val;
        double ln_1m_p_ir = std::log(1-(**p_ir));
        int Rstar_val, Sstar_val, Estar_val, R_val, I_val, S_val;   

        if (((*A0) -> I0)[startLoc] < 0 || 
            ((*A0) -> R0)[startLoc] < 0)
        {
            lssCout << "I0 < 0: " << (((*A0) -> I0)[startLoc]) << "\n"; 
            throw(-1);
        }
        compIdx = startLoc*nTpts;

        // Is p_rs meaningful?
        if ((*context) -> config -> reinfectionMode <= 2)
        {

            for (j = 0; j < nTpts; j++)
            {
                Rstar_val = ((*R_star) -> data)[compIdx];
                Sstar_val = ((*S_star)->data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];
                p_rs_val = (*p_rs)[j];

                if (Rstar_val > I_val || 
                        Sstar_val > R_val)
                {
                    lssCout << "Bounds error detected at time " << j << "\n";
                    lssCout << "R_star: " << Rstar_val << "\n";
                    lssCout << "I: " << I_val << "\n";
                    lssCout << "S_star: " << Sstar_val << "\n";
                    lssCout << "R: " << R_val << "\n";
                    throw(-1) ;
                }
                else
                {
                    output +=  (ln_1m_p_ir*(I_val) +
                                ((*context) -> random -> dbinom(Sstar_val, R_val, p_rs_val)) + 
                                ((*context) -> random -> choose(I_val, Rstar_val)));
                    if (!std::isfinite(output))
                    {
                        lssCout << "Computation error detected at time " << j << "\n";
                        lssCout << "S_star: " << Sstar_val << "\n";
                        lssCout << "R_star: " << Rstar_val << "\n";
                        lssCout << "R: " << R_val << "\n";
                        lssCout << "I: " << I_val << "\n";
                        lssCout << "p_rs: " << p_rs_val << "\n";
                        throw(-1);
                    }
                }
                compIdx++;
            } 
        }
        else
        {
            for (j = 0; j < nTpts; j++)
            {
                Rstar_val = ((*R_star) -> data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];

                if (Rstar_val > I_val)
                {
                    lssCout << "Bounds error detected at time " << j << "\n";
                    throw(-1);
                }
                else
                {
                    output +=  (ln_1m_p_ir*(I_val) +
                                ((*context) -> random -> choose(I_val, Rstar_val)));
                    if (!std::isfinite(output))
                    {
                        lssCout << "Computation error detected at time " << j << "\n";
                        lssCout << "R_star: " << Rstar_val << "\n";
                        lssCout << "R: " << R_val << "\n";
                        lssCout << "I: " << I_val << "\n";
                        throw(-1);
                    }

                }
                compIdx++;
            } 

        }

        // p_se changes, so need to look at p_se component for all locations and 
        // time points after 0
        for (i = 0; i < nLoc; i++)
        {
            compIdx = i*nTpts;
            for (j = 0; j< nTpts; j++)
            {

                p_se_val = (*p_se)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                if (p_se_val > 1 || p_se_val < 0)
                {
                    lssCout << "p_se bounds error detected at time " << j << ": " << p_se_val << "\n";
                    throw(-1);
                }

                output += (*context) -> random -> dbinom(Estar_val,S_val, p_se_val);
                if (!std::isfinite(output))
                {
                    lssCout << "p_se computation error detected at time " << j << "\n";
                    lssCout << "S: " << S_val << "\n";
                    lssCout << "E_star: " << Estar_val << "\n";
                    lssCout << "p_se: " << p_se_val << "\n";
                }
                compIdx ++; 
            }
        }

        if (!std::isfinite(output))
        {

            throw(-1);
        }
        return;
    }

}
