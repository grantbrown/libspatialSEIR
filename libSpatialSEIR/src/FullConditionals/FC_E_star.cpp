#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_FC_E_star.hpp>
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
                         double _sliceWidth,
                         int _useOCL) 
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
        useOCL = new int;
       
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
        *useOCL = _useOCL;
    }

    FC_E_Star::~FC_E_Star()
    {
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
        delete useOCL;

    }

    int FC_E_Star::evalCPU(int startLoc, int startTime)
    {
        int j, compIdx;
        int nTpts = *((*S) -> nrow); 
        double ln_1m_p_ei = std::log(1-(**p_ei));    
        double p_se_val;
        int S_val, E_val, Estar_val, Istar_val;
        long double output = 0.0;
        long unsigned int E_star_sum;
        long unsigned int I_star_sum;
        int64_t aDiff; 

        compIdx = startLoc*nTpts + startTime;
        for (j = startTime; j < nTpts; j++)
        {
            Estar_val = ((*E_star) -> data)[compIdx];
            S_val = ((*S) -> data)[compIdx];
            E_val = ((*E) -> data)[compIdx];
            Istar_val = ((*I_star) -> data)[compIdx];
            p_se_val = (*p_se)[compIdx];

            if (Estar_val < 0 || Estar_val > S_val || 
                    Istar_val > E_val)

            {
                *value = -INFINITY;
                return(-1);
            }
            else
            {
                output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)) +    
                             ln_1m_p_ei*E_val + 
                             ((*context) -> random -> choosePartial(E_val, Istar_val)));
            }
            compIdx++;
        }

        E_star_sum = (*E_star)->marginSum(2,startLoc);
        I_star_sum = (*I_star)->marginSum(2,startLoc);
        aDiff = (E_star_sum > I_star_sum ? E_star_sum - I_star_sum : I_star_sum - E_star_sum)/nTpts;
        output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);

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

    int FC_E_Star::evalCPU()
    {
        int i,j, compIdx;
        int nTpts = *((*S) -> nrow); 
        int nLoc = *((*S) -> ncol); 
        double ln_1m_p_ei = std::log(1-(**p_ei));    
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

                if (Estar_val < 0 || Estar_val > S_val || 
                        Istar_val > E_val)

                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                {
                    output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)) +    
                                 ln_1m_p_ei*E_val + 
                                 ((*context) -> random -> choosePartial(E_val, Istar_val)));
                }
                compIdx++;
            }

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

    void FC_E_Star::printDebugInfo(int loc, int tpt)
    {
        std::cout << "E_star debug info, location " << loc << ", time " << tpt << "\n";     
        std::cout << "...looking for problems\n";

        int j, compIdx;
        int nTpts = *((*S) -> nrow); 
        double ln_1m_p_ei = std::log(1-(**p_ei));    
        double p_se_val;
        int S_val, E_val, Estar_val, Istar_val;
        long double output = 0.0;
        long unsigned int E_star_sum;
        long unsigned int I_star_sum;
        int64_t aDiff; 


        compIdx = loc*nTpts + tpt;
        for (j = tpt; j < nTpts; j++)
        {
            std::cout << "time " << j << "\n";
            Estar_val = ((*E_star) -> data)[compIdx];
            S_val = ((*S) -> data)[compIdx];
            E_val = ((*E) -> data)[compIdx];
            Istar_val = ((*I_star) -> data)[compIdx];
            p_se_val = (*p_se)[compIdx];

            if (Estar_val < 0 || Estar_val > S_val || 
                    Istar_val > E_val)

            {
                std::cout << "Bounds Error Detected:\n";
                std::cout << "E_star: " << Estar_val << "\n";
                std::cout << "S: " << S_val << "\n";
                std::cout << "I_star: " << Istar_val << "\n";
                return;
            }
            else
            {
                output = (std::log(p_se_val)*Estar_val +
                             std::log(1-p_se_val)*(S_val - Estar_val) +
                             ln_1m_p_ei*E_val + 
                             ((*context) -> random -> choosePartial(E_val, Istar_val)) + 
                             ((*context)->random->choose(S_val, Estar_val)));
                if (!std::isfinite(output))
                {
                    std::cout << "Computation Error Detected: " << tpt << "\n";
                    std::cout << "E_star: " << Estar_val << "\n";
                    std::cout << "S: " << S_val << "\n";
                    std::cout << "E: " << E_val << "\n";
                    std::cout << "I_star: " << Istar_val << "\n";
                    std::cout << "p_se: " << p_se_val << "\n";
                    std::cout << "ln_1m_p_ei: " << ln_1m_p_ei << "\n";
                    return; 
                }
            }
            compIdx++;
        }

        E_star_sum = (*E_star)->marginSum(2,loc);
        I_star_sum = (*I_star)->marginSum(2,loc);
        aDiff = (E_star_sum > I_star_sum ? E_star_sum - I_star_sum : I_star_sum - E_star_sum)/nTpts;
        output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);

        if (!std::isfinite(output))
        {
            std::cout << "Combinatorics Error Detectd:\n";
            std::cout << "Computation Error Detected:\n";
            std::cout << "E_star: " << Estar_val << "\n";
            std::cout << "S: " << S_val << "\n";
            std::cout << "E: " << E_val << "\n";
            std::cout << "I_star: " << Istar_val << "\n";
            std::cout << "p_se: " << p_se_val << "\n";
            std::cout << "ln_1m_p_ei: " << ln_1m_p_ei << "\n";
            std::cout << "E_star_sum: " << E_star_sum << "\n";
            std::cout << "I_star_sum: " << I_star_sum << "\n";
            return;
        }    
       return;
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
    int FC_E_Star::calculateRelevantCompartments_OCL()
    {
        (*context) -> calculateE_CPU();
        (*context) -> calculateS_givenE_CPU();
        return(0);
    }
    int FC_E_Star::calculateRelevantCompartments(int startLoc, int startTime)
    {
        (*context) -> calculateE_CPU(startLoc, startTime);
        (*context) -> calculateS_givenE_CPU(startLoc, startTime);
        return(0);
    }

    void FC_E_Star::sample(int verbose)
    {
        if (verbose){std::cout << "Sampling E_star\n";}
        if (*useOCL){sampleOCL(); return;}
        sampleCPU();
    }

    int FC_E_Star::sampleCPU()
    {
        int mode = (*context) -> getSamplingMode();
        if (mode == 1)
        {
            this -> sampleCompartment_CPU(*context,
                                          *E_star,*sliceWidth);
        }
        else if (mode == 2)
        {
            this -> sampleEntireCompartment_CPU(*context,
                                       *E_star,*sliceWidth);
        }
        else if (mode == 3)
        {
            this -> sampleEntireCompartment2_CPU(*context,
                                       *E_star,*sliceWidth);

        }
        else if (mode == 4)
        {
            if ((*context) -> random -> uniform() < 0.9)
            {
                this -> sampleEntireCompartment2_CPU(*context,
                                          *E_star,*sliceWidth);
            }
            else 
            {
                this -> sampleCompartment_CPU(*context,
                                          *E_star,*sliceWidth);
            }
        }
        else
        {
            std::cout << "Invalid sampling mode, falling back to default.\n";
            this -> sampleEntireCompartment_CPU(*context,
                                       *E_star,*sliceWidth);
        }

        return 0;
    }

    int FC_E_Star::sampleOCL()
    {
        this -> sampleCompartment_OCL(*context,
                                  *E_star,*sliceWidth);
        return(0);
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
