#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FC_S_star.hpp>
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
     * Implement the full conditional distribution for S_star
     *
     */    



    FC_S_Star::FC_S_Star(ModelContext * _context,
                         CompartmentalModelMatrix *_S_star, 
                         CompartmentalModelMatrix *_S, 
                         CompartmentalModelMatrix *_R, 
                         CompartmentalModelMatrix *_E_star,
                         CompartmentalModelMatrix *_R_star,
                         InitData *_A0,
                         CovariateMatrix *_X, 
                         double *_p_se,
                         double *_p_rs,
                         double *_beta,
                         double *_rho,
                         double _steadyStateConstraintPrecision,
                         double _sliceWidth)
    {
       context = new ModelContext*;
       S_star = new CompartmentalModelMatrix*;
       S = new CompartmentalModelMatrix*;
       R = new CompartmentalModelMatrix*;
       E_star = new CompartmentalModelMatrix*;
       R_star = new CompartmentalModelMatrix*;
       A0 = new InitData*;
       X = new CovariateMatrix*;
       p_se = new double*;
       p_rs = new double*;
       beta = new double*;
       rho = new double*;
       value = new long double;
       steadyStateConstraintPrecision = new double;
       sliceWidth = new double;
       samples = new int;
       accepted = new int; 
       *samples = 0;
       *accepted = 0;

       *context = _context;
       *S_star = _S_star;
       *S = _S;
       *R = _R;
       *E_star = _E_star;
       *R_star = _R_star;
       *A0 = _A0;
       *X = _X;
       *p_se = _p_se;
       *p_rs = _p_rs;
       *beta = _beta;
       *rho = _rho;
       *steadyStateConstraintPrecision = _steadyStateConstraintPrecision;
       *value = -1.0;
       *sliceWidth = _sliceWidth;

        // Set up samplers
        samplers = new std::vector<Sampler*>();
        currentSampler = new Sampler*;
        samplers -> push_back(new CompartmentMetropolisSampler(*context, this, (*S_star) -> data));
        samplers -> push_back(new IndexedCompartmentMetropolisSampler(*context, this, (*S_star) -> data));
        samplers -> push_back(new CompartmentMetropolisSampler_OCL(*context, this, (*S_star) -> data));


    }    
    FC_S_Star::~FC_S_Star()
    {
        while((samplers -> size()) != 0){delete (*samplers).back(); (*samplers).pop_back();}
        delete samplers;
        delete S_star;
        delete S;
        delete R;
        delete E_star;
        delete R_star;
        delete A0;
        delete X;
        delete p_se;
        delete p_rs;
        delete beta;
        delete rho;
        delete value;
        delete steadyStateConstraintPrecision;
        delete sliceWidth;
        delete context;
        delete samples;
        delete accepted;
    }

    int FC_S_Star::evalCPU(int startLoc, int startTime)
    {
        int i,compIdx,Sstar_val,Estar_val,S_val,R_val;
        double p_se_val, p_rs_val;
        int nTpts = *((*S)->nrow);
        long double output = 0.0;
        long unsigned int S_star_sum;
        long unsigned int R_star_sum;
        int64_t aDiff;

        compIdx = startLoc*nTpts + startTime;
        for (i = startTime; i < nTpts; i++)
        {
                Sstar_val = ((*S_star)->data)[compIdx]; 
                Estar_val = ((*E_star)->data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                R_val = ((*R)->data)[compIdx];
                p_se_val = (*p_se)[compIdx];
                p_rs_val = (*p_rs)[i];

                if (Sstar_val < 0 || 
                    Sstar_val > R_val ||
                    Estar_val > S_val)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                { 
                    output += (((*context) -> random -> dbinom(Sstar_val, R_val, p_rs_val)) + 
                               ((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)));
                }
                compIdx ++; 
        }

        if (*steadyStateConstraintPrecision > 0)
        {
            S_star_sum = (*S_star)->marginSum(2,startLoc);
            R_star_sum = (*R_star)->marginSum(2,startLoc);
            aDiff = (S_star_sum > R_star_sum ? S_star_sum - R_star_sum : R_star_sum - S_star_sum)/nTpts;
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
        return(0);
    }

    int FC_S_Star::evalCPU()
    {
        int i,j,compIdx,Sstar_val,Estar_val,S_val,R_val;
        double p_se_val, p_rs_val;
        int nTpts = *((*S)->nrow);
        int nLoc = *((*S)->ncol);
        long double output = 0.0;
        long unsigned int S_star_sum;
        long unsigned int R_star_sum;
        int64_t aDiff;

        for (i = 0; i < nLoc ; i++)
        {
            compIdx = i*nTpts;
            
            for (j = 0; j < nTpts; j++)
            {
                    Sstar_val = ((*S_star)->data)[compIdx]; 
                    Estar_val = ((*E_star)->data)[compIdx];
                    S_val = ((*S)->data)[compIdx];
                    R_val = ((*R)->data)[compIdx];
                    p_se_val = (*p_se)[compIdx];
                    p_rs_val = (*p_rs)[j];

                    if (Sstar_val < 0 || 
                        Sstar_val > R_val ||
                        Estar_val > S_val)
                    {
                        *value = -INFINITY;
                        return(-1);
                    }
                    else
                    { 
                        output += (((*context) -> random -> dbinom(Sstar_val, R_val, p_rs_val)) + 
                                   ((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)));
                    }
                    compIdx ++; 
            }

            if (*steadyStateConstraintPrecision > 0)
            {
                S_star_sum = (*S_star)->marginSum(2,i);
                R_star_sum = (*R_star)->marginSum(2,i);
                aDiff = (S_star_sum > R_star_sum ? S_star_sum - R_star_sum : R_star_sum - S_star_sum)/nTpts;
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
        return(0);
    }

    void FC_S_Star::printDebugInfo(int startLoc, int startTime)
    {
        lssCout << "S_star debug info, location " << startLoc << ", time " << startTime << "\n";     

        int i,compIdx,Sstar_val,Estar_val,S_val,R_val;
        double p_se_val, p_rs_val;
        int nTpts = *((*S)->nrow);
        long double output = 0.0;
        long unsigned int S_star_sum;
        long unsigned int R_star_sum;
        int64_t aDiff;

        compIdx = startLoc*nTpts + startTime;
        for (i = startTime; i < nTpts; i++)
        {
                Sstar_val = ((*S_star)->data)[compIdx]; 
                Estar_val = ((*E_star)->data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                R_val = ((*R)->data)[compIdx];
                p_se_val = (*p_se)[compIdx];
                p_rs_val = (*p_rs)[i];

                if (Sstar_val < 0 || 
                    Sstar_val > R_val ||
                    Estar_val > S_val)
                {
                    lssCout << "Bounds Error Detected at time " << i << "\n";
                    lssCout << "S_star: " << Sstar_val << "\n";
                    lssCout << "E_star: " << Estar_val << "\n";
                    lssCout << "R: " << R_val << "\n";
                    lssCout << "S: " << S_val << "\n";
                    return;
                }
                else
                { 
                    output += (((*context) -> random -> dbinom(Sstar_val, R_val, p_rs_val)) + 
                               ((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)));
                }
                if (! std::isfinite(output))
                {
                    lssCout << "Calculation Error Detected at time " << i << "\n";
                    lssCout << "S_star: " << Sstar_val << "\n";
                    lssCout << "E_star: " << Estar_val << "\n";
                    lssCout << "R: " << R_val << "\n";
                    lssCout << "S: " << S_val << "\n";
                    lssCout << "p_se: " << p_se_val << "\n";
                    lssCout << "p_rs: " << p_rs_val << "\n"; 
                    return;
                }
                compIdx ++; 
        }
        if (*steadyStateConstraintPrecision > 0)
        {
            S_star_sum = (*S_star)->marginSum(2,startLoc);
            R_star_sum = (*R_star)->marginSum(2,startLoc);
            aDiff = (S_star_sum > R_star_sum ? S_star_sum - R_star_sum : R_star_sum - S_star_sum)/nTpts;
            output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);
        }

        if (!std::isfinite(output))
        {
            lssCout << "Combinatorics Error Detected\n";
            lssCout << "S_star: " << Sstar_val << "\n";
            lssCout << "E_star: " << Estar_val << "\n";
            lssCout << "R: " << R_val << "\n";
            lssCout << "S: " << S_val << "\n";
            lssCout << "p_se: " << p_se_val << "\n";
            lssCout << "p_rs: " << p_rs_val << "\n"; 
            lssCout << "S_star_sum: " << S_star_sum << "\n";
            lssCout << "R_star_sum: " << R_star_sum << "\n";
            return;
        }
        return;

    }


    int FC_S_Star::evalOCL()
    {
        // Not Implemented
        return(evalCPU());
    }
    int FC_S_Star::calculateRelevantCompartments()
    {
        (*context) -> calculateS_CPU();
        (*context) -> calculateR_givenS_CPU();
        return(0);
    }
    int FC_S_Star::calculateRelevantCompartments_OCL()
    {
        (*context) -> calculateS_CPU();
        (*context) -> calculateR_givenS_CPU();
        return(0);

    }
    int FC_S_Star::calculateRelevantCompartments(int startLoc, int startTime)
    {
        (*context) -> calculateS_CPU(startLoc, startTime);
        (*context) -> calculateR_givenS_CPU(startLoc, startTime);
        return(0);
    }

    void FC_S_Star::sample(int verbose)
    {
        if (verbose){lssCout << "Sampling S_star\n";}
        (*currentSampler) -> drawSample();
    }

    long double FC_S_Star::getValue()
    {
        return(*(this -> value));
    }
    void FC_S_Star::setValue(long double val)
    {
        *(this -> value) = val;
    }
}
