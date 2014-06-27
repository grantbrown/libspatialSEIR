#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_FC_E0.hpp>
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
     * Implement full conditional for E0
     *
     */
    FC_E0::FC_E0(ModelContext* _context, 
                 CompartmentalModelMatrix *_S,
                 CompartmentalModelMatrix *_E,
                 CompartmentalModelMatrix *_I,
                 CompartmentalModelMatrix *_E_star,
                 CompartmentalModelMatrix *_I_star,
                 CompartmentalModelMatrix *_R_star,
                 InitData *_A0,
                 double *_p_ir,
                 double *_p_ei,
                 double *_p_se,
                 double _sliceWidth)
    {
        context = new ModelContext*;
        S = new CompartmentalModelMatrix*;
        E = new CompartmentalModelMatrix*;
        I = new CompartmentalModelMatrix*;
        E_star = new CompartmentalModelMatrix*;
        I_star = new CompartmentalModelMatrix*;
        R_star = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_se = new double*;
        p_ir = new double*;
        p_ei = new double*;
        sliceWidth = new double;
        value = new long double;
        samples = new int; *samples = 0;
        accepted = new int; *accepted =0; 


        *context = _context;
        *S = _S;
        *E = _E;
        *I = _I;
        *E_star = _E_star;
        *I_star = _I_star;
        *R_star = _R_star;
        *p_se = _p_se;
        *p_ir = _p_ir;
        *p_ei = _p_ei;
        *A0 = _A0;
        *sliceWidth = _sliceWidth;
    }
    FC_E0::~FC_E0()
    {
        delete context;
        delete S;
        delete E;
        delete I;
        delete E_star;
        delete I_star;
        delete R_star;
        delete p_se;
        delete p_ir;
        delete p_ei;
        delete A0;
        delete sliceWidth;
        delete value;
        delete samples;
        delete accepted;

    }

    int FC_E0::evalCPU()
    {
        int i,j,compIdx,S_val,E_val,I_val,Istar_val,Estar_val,Rstar_val;
        double p_ei_val, p_ir_val, p_se_val;
        int nLoc = *((*E)->ncol);
        int nTpts = *((*E)->nrow);
        long double output = 0.0;

        p_ei_val = **p_ei;
        p_ir_val = **p_ir;
        for (i = 0; i < nLoc; i++)
        {
            if (((*A0) -> E0)[i] < 0 || 
                ((*A0) -> I0)[i] < 0)
            {
                *value = -INFINITY;
                return(-1);
            }

            compIdx = i*nTpts;
            for (j = 0; j < nTpts; j++)
            {
                Rstar_val = ((*R_star)->data)[compIdx]; 
                Estar_val = ((*E_star) -> data)[compIdx];
                Istar_val = ((*I_star)->data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                E_val = ((*E)->data)[compIdx];
                I_val = ((*I)->data)[compIdx];
                p_se_val = (*p_se)[compIdx];

                if (Istar_val > E_val ||
                    Rstar_val > I_val || 
                    p_se_val > 1 || 
                    p_se_val < 0)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                { 
                    output += (((*context) -> random -> dbinom(Rstar_val, I_val, p_ir_val)) + 
                               ((*context) -> random -> dbinom(Istar_val, E_val, p_ei_val)) + 
                               ((*context) -> random -> dbinom(Estar_val,S_val, p_se_val)));

                }
                compIdx ++; 
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
        return(0);
    }

    
    int FC_E0::evalCPU(int startLoc)
    {
        int i,j,compIdx,S_val,E_val,I_val,Istar_val,Estar_val,Rstar_val;
        double p_ei_val, p_ir_val, p_se_val;
        int nLoc = *((*E)->ncol);
        int nTpts = *((*E)->nrow);
        long double output = 0.0;

        if (((*A0) -> E0)[startLoc] < 0 || 
            ((*A0) -> I0)[startLoc] < 0)
        {
            *value = -INFINITY;
            return(-1);
        }

        compIdx = startLoc*nTpts;
        p_ei_val = **p_ei;
        p_ir_val = **p_ir;
        for (i = 0; i < nTpts; i++)
        {
                Rstar_val = ((*R_star)->data)[compIdx]; 
                Istar_val = ((*I_star)->data)[compIdx];
                E_val = ((*E)->data)[compIdx];
                I_val = ((*I)->data)[compIdx];
                if (Istar_val > E_val ||
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
                compIdx ++; 
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
                    return(-1);
                }

                output += (*context) -> random -> dbinom(Estar_val,S_val, p_se_val);
                compIdx ++; 
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
        return(0);
    }

    int FC_E0::evalOCL()
    {
        // Not Implemented
        return(evalCPU());
    }
    int FC_E0::sampleCPU()
    {
        sampleEntireCompartment_CPU(*context, (*A0) -> E0, *sliceWidth);
        return(0);
    }
    int FC_E0::sampleOCL()
    {
        return(sampleEntireCompartment_OCL(*context, (*A0) -> E0, *sliceWidth));
    }

    long double FC_E0::getValue()
    {
        return(*value);
    }

    void FC_E0::setValue(long double val)
    {
        *(this -> value) = val;
    }

    int FC_E0::calculateRelevantCompartments()
    {
        (*context) -> calculateE_CPU();
        (*context) -> calculateI_givenR_CPU();
        (*context) -> calculateP_SE_CPU();
        return(0);
    }

    int FC_E0::calculateRelevantCompartments_OCL()
    {
        (*context) -> calculateE_CPU();
        (*context) -> calculateI_givenR_CPU();
        (*context) -> calculateP_SE_OCL();
        return(0);
    }

    int FC_E0::calculateRelevantCompartments(int startLoc)
    {
        (*context) -> calculateE_CPU(startLoc, 0);
        (*context) -> calculateI_givenR_CPU(startLoc,0);
        (*context) -> calculateP_SE_CPU(startLoc,0);
        return(0);
    }

    void FC_E0::printDebugInfo(int loc)
    {
        std::cout << "Error Sampling E0, location: " << loc << ", value: " << ((*A0) -> E0)[loc] << "\n";
    }
}
