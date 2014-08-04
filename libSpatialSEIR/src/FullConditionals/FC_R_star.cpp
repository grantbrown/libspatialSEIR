#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FC_R_star.hpp>
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
     * Implement the full conditional distribution for R_Star
     *
     */    
    FC_R_Star::FC_R_Star(ModelContext *_context,
                         CompartmentalModelMatrix *_R_star,
                         CompartmentalModelMatrix *_R,
                         CompartmentalModelMatrix *_I,
                         CompartmentalModelMatrix *_S_star,
                         CompartmentalModelMatrix *_E_star,
                         CompartmentalModelMatrix *_I_star,
                         CompartmentalModelMatrix *_S,
                         InitData *_A0,
                         double *_p_rs,
                         double *_p_ir,
                         double *_p_se,
                         double _steadyStateConstraintPrecision,
                         double _sliceWidth)
    {

        context = new ModelContext*;
        R_star = new CompartmentalModelMatrix*;
        R = new CompartmentalModelMatrix*;
        I = new CompartmentalModelMatrix*;
        S_star = new CompartmentalModelMatrix*;
        E_star = new CompartmentalModelMatrix*;
        I_star = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_rs = new double*;
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
        *R_star = _R_star;
        *R = _R;
        *I = _I;
        *S_star = _S_star;
        *E_star = _E_star;
        *I_star = _I_star;
        *S = _S;
        *A0 = _A0;
        *p_rs = _p_rs;
        *p_ir = _p_ir;
        *p_se = _p_se;
        *sliceWidth = _sliceWidth;
        *steadyStateConstraintPrecision = _steadyStateConstraintPrecision;
        *value = -1.0;

        // Set up samplers
        samplers = new std::vector<Sampler*>();
        currentSampler = new Sampler*;
        samplers -> push_back(new CompartmentMetropolisSampler(*context, this, (*R_star) -> data));
        samplers -> push_back(new IndexedCompartmentMetropolisSampler(*context, this, (*R_star) -> data));
        samplers -> push_back(new CompartmentMetropolisSampler_OCL(*context, this, (*R_star) -> data));

    }
    FC_R_Star::~FC_R_Star()
    {
        while((samplers -> size()) != 0){delete (*samplers).back(); (*samplers).pop_back();}

        delete samplers;
        delete R_star;
        delete R;
        delete I;
        delete S_star;
        delete E_star;
        delete I_star;
        delete S;
        delete A0;
        delete p_rs;
        delete p_ir;
        delete p_se;
        delete value;
        delete sliceWidth;
        delete steadyStateConstraintPrecision;
        delete context;
        delete samples;
        delete accepted;
    }

    int FC_R_Star::evalCPU(int startLoc, int startTime)
    {
        int i,j, compIdx;
        int nTpts = *((*R) -> nrow);
        int nLoc = *((*R) -> ncol);

        long double output = 0.0;
        
        double p_se_val;
        double p_rs_val;
        double p_ir_val;
        int Rstar_val, Sstar_val, Estar_val, R_val, I_val, S_val;   
        long unsigned int I_star_sum;
        long unsigned int R_star_sum;
        int64_t aDiff; 

        compIdx = startLoc*nTpts + startTime;
        // Is p_rs meaningful?
        if ((*context) -> config -> reinfectionMode <= 2)
        {
            for (j = startTime; j < nTpts; j++)
            {
                Rstar_val = ((*R_star) -> data)[compIdx];
                Sstar_val = ((*S_star)->data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];
                p_rs_val = (*p_rs)[j];
                p_ir_val = (*p_ir)[j];
                if (Rstar_val < 0 || Rstar_val > I_val || 
                        Sstar_val > R_val)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                {
                    output += (((*context) -> random -> dbinom(Rstar_val, I_val, p_ir_val)) + 
                               ((*context) -> random -> dbinom(Sstar_val, R_val, p_rs_val)));
                }
                compIdx++;
            } 
        }
        else
        {
            for (j = startTime; j < nTpts; j++)
            {
                p_ir_val = (*p_ir)[j];
                Rstar_val = ((*R_star) -> data)[compIdx];
                I_val = ((*I) ->data)[compIdx];

                if (Rstar_val < 0 || Rstar_val > I_val)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                {
                    output += ((*context) -> random -> dbinom(Rstar_val, I_val, p_ir_val)); 
                }
                compIdx++;
            } 

        }

        // p_se changes, so need to look at p_se component for all locations and 
        // time points after startTime
        for (i = 0; i < nLoc; i++)
        {
            compIdx = i*nTpts + startTime;
            for (j = startTime; j< nTpts; j++)
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

        if (*steadyStateConstraintPrecision > 0)
        {
            I_star_sum = (*I_star)->marginSum(2,startLoc);
            R_star_sum = (*R_star)->marginSum(2,startLoc);
            aDiff = (I_star_sum > R_star_sum ? I_star_sum - R_star_sum : R_star_sum - I_star_sum)/nTpts;
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

    int FC_R_Star::evalCPU()
    {
        int i,j, compIdx;
        int nTpts = *((*R) -> nrow);
        int nLoc = *((*R) -> ncol);

        long double output = 0.0;
        
        double p_se_val;
        double p_rs_val;
        double p_ir_val;
        int Rstar_val, Sstar_val, Estar_val, R_val, I_val, S_val;   
        long unsigned int I_star_sum;
        long unsigned int R_star_sum;
        int64_t aDiff; 

        // Is p_rs meaningful?
        if ((*context) -> config -> reinfectionMode <= 2)
        {
            for (i = 0; i < nLoc; i++)
            {
                compIdx = i*nTpts;
                for (j = 0; j < nTpts; j++)
                {
                    Rstar_val = ((*R_star) -> data)[compIdx];
                    Sstar_val = ((*S_star)->data)[compIdx];
                    R_val = ((*R) ->data)[compIdx];
                    I_val = ((*I) ->data)[compIdx];
                    p_rs_val = (*p_rs)[j];
                    p_ir_val = (*p_ir)[j];

                    if (Rstar_val < 0 || Rstar_val > I_val || 
                            Sstar_val > R_val)
                    {
                        *value = -INFINITY;
                        return(-1);
                    }
                    else
                    {
                        output += (((*context) -> random -> dbinom(Rstar_val, I_val, p_ir_val)) + 
                                   ((*context) -> random -> dbinom(Sstar_val, R_val, p_rs_val)));
                    }
                    compIdx++;
                } 
            }
        }
        else
        {
            for (i = 0; i < nLoc; i++)
            {
                compIdx = i*nTpts;
                for (j = 0; j < nTpts; j++)
                {
                    Rstar_val = ((*R_star) -> data)[compIdx];
                    I_val = ((*I) ->data)[compIdx];
                    p_ir_val = (*p_ir)[j];

                    if (Rstar_val < 0 || Rstar_val > I_val)
                    {
                        *value = -INFINITY;
                        return(-1);
                    }
                    else
                    {
                        output += ((*context) -> random -> dbinom(Rstar_val, I_val, p_ir_val));
                    }
                    compIdx++;
                } 
            }
        }


        // p_se changes, so need to look at p_se component for all locations and 
        // time points
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

        if (*steadyStateConstraintPrecision > 0)
        {
            I_star_sum = (*I_star)->marginSum(3,-1);
            R_star_sum = (*R_star)->marginSum(3,-1);
            aDiff = (I_star_sum > R_star_sum ? I_star_sum - R_star_sum : R_star_sum - I_star_sum)/(nTpts*nLoc);
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


    void FC_R_Star::printDebugInfo(int startLoc, int startTime)
    {
        lssCout << "R_star debug info, location " << startLoc << ", time " << startTime << "\n";     
        lssCout << "...looking for problems\n";


        int i,j, compIdx;
        int nTpts = *((*R) -> nrow);
        int nLoc = *((*R) -> ncol);

        long double output = 0.0;
        
        double p_se_val;
        double p_rs_val;
        double p_ir_val;
        int Rstar_val, Sstar_val, Estar_val, R_val, I_val, S_val;   
        long unsigned int I_star_sum;
        long unsigned int R_star_sum;
        int64_t aDiff; 

        compIdx = startLoc*nTpts + startTime;
        if ((*context) -> config -> reinfectionMode <= 2)
        {
            for (j = startTime; j < nTpts; j++)
            {
                Rstar_val = ((*R_star) -> data)[compIdx];
                Sstar_val = ((*S_star)->data)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                p_rs_val = (*p_rs)[j];
                p_ir_val = (*p_ir)[j];

                if (Rstar_val < 0 || Rstar_val > I_val || 
                        Sstar_val > R_val)
                {

                    lssCout << "Bounds error detected, time " << j << "\n";
                    lssCout << "S_star: " << Sstar_val << "\n";
                    lssCout << "R_star: " << Rstar_val << "\n";
                    lssCout << "I: " << I_val << "\n";
                    lssCout << "R: " << R_val << "\n";
                    return;
                }
                else
                {
                    output += ( ((*context) -> random -> dbinom(Sstar_val, R_val, p_rs_val)) + 
                                ((*context) -> random -> dbinom(Rstar_val, I_val, p_ir_val)));
                    if (!std::isfinite(output))
                    {
                        lssCout << "Calculation Error Detected, time " << j << "\n";
                        lssCout << "S_star: " << Sstar_val << "\n";
                        lssCout << "R_star: " << Rstar_val << "\n";
                        lssCout << "I: " << I_val << "\n";
                        lssCout << "R: " << R_val << "\n";
                        lssCout << "p_ir: " << p_ir_val << "\n";
                        return;
                    }
                }
                compIdx++;
            } 
        }
        else 
        {
            for (j = startTime; j < nTpts; j++)
            {
                Rstar_val = ((*R_star) -> data)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                p_ir_val = (*p_ir)[j];

                if (Rstar_val < 0 || Rstar_val > I_val)
                {

                    lssCout << "Bounds error detected, time " << j << "\n";
                    lssCout << "R_star: " << Rstar_val << "\n";
                    lssCout << "I: " << I_val << "\n";
                    lssCout << "R: " << R_val << "\n";
                    return;
                }
                else
                {
                    output += ((*context) -> random -> dbinom(Rstar_val, I_val, p_ir_val));
                    if (!std::isfinite(output))
                    {
                        lssCout << "Calculation Error Detected, time " << j << "\n";
                        lssCout << "R_star: " << Rstar_val << "\n";
                        lssCout << "I: " << I_val << "\n";
                        lssCout << "R: " << R_val << "\n";
                        lssCout << "p_ir: " << p_ir_val << "\n";
                        return;
                    }
                }
                compIdx++;
            } 

        }

        // p_se changes, so need to look at p_se component for all locations and 
        // time points after startTime
        for (i = 0; i < nLoc; i++)
        {
            compIdx = i*nTpts + startTime;
            for (j = startTime; j< nTpts; j++)
            {
                p_se_val = (*p_se)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                if (p_se_val > 1 || p_se_val < 0)
                {
                    lssCout << "Non-Local Bounds Error Detected, loc, time: " << i << ", " << j << "\n";
                    lssCout << "E_star: " << Estar_val << "\n";
                    lssCout << "S: " << S_val << "\n";
                    lssCout << "p_se: " << p_se_val << "\n";
                    return;
                }

                output += (*context) -> random -> dbinom(Estar_val, S_val, p_se_val);

                if (!std::isfinite(output))
                {
                    lssCout<< "Non-Local Calculation Error Detected, loc, time: " << i << ", " << j << "\n";
                    lssCout << "E_star: " << Estar_val << "\n";
                    lssCout << "S: " << S_val << "\n";
                    lssCout << "p_se: " << p_se_val << "\n";
                    lssCout << "log(p_se): " << std::log(p_se_val) << "\n";
                    lssCout << "p_se == 0: " << (p_se_val == 0.0) << "\n";
                    lssCout << "Binomial term: " << ((*context) -> random -> 
                            dbinom(Estar_val, S_val, p_se_val)) << "\n";


                    return;
                }

                compIdx ++; 
            }
        }

        if (*steadyStateConstraintPrecision > 0)
        {
            I_star_sum = (*I_star)->marginSum(2,startLoc);
            R_star_sum = (*R_star)->marginSum(2,startLoc);
            aDiff = (I_star_sum > R_star_sum ? I_star_sum - R_star_sum : R_star_sum - I_star_sum)/nTpts;
            output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);
        }

        if (!std::isfinite(output))
        {
            lssCout << "Combinatorics Error Detected.\n";
            return;
        }

        lssCout << "No problems detected?\n";
        return;
    }


    int FC_R_Star::evalOCL()
    {

        int nTpts = *((*R) -> nrow);
        int nLoc = *((*R) -> ncol);

        if ((*context) -> config -> reinfectionMode > 2)
        {
            lssCout << "FC_R_Star currently only works with OpenCL for reinfectionMode <= 2\n";
        }

        double output = ((*context) -> oclProvider -> 
                FC_R_Star(nLoc,
                          nTpts,
                          ((*S_star) -> data),
                          ((*E_star) -> data),
                          ((*R_star) -> data),
                          ((*S) -> data),
                          ((*I) -> data),
                          ((*R) -> data),
                          (*p_se),
                          (*p_rs),
                          (*p_ir)
                          ));
        if (!std::isfinite(output))
        {
            *value = -INFINITY;
            return(-1);
        }
        else 
        *value = output;
        return 0;
    }
    int FC_R_Star::calculateRelevantCompartments()
    {        
        (*context) -> calculateR_CPU();
        (*context) -> calculateI_givenR_CPU();
        ((*context) -> calculateP_SE_CPU());
        return(0);
    }

    int FC_R_Star::calculateRelevantCompartments_OCL()
    {
        (*context) -> calculateR_CPU();
        (*context) -> calculateI_givenR_CPU();
        ((*context) -> calculateP_SE_OCL());
        return(0);
    }

    int FC_R_Star::calculateRelevantCompartments(int startLoc, int startTime)
    {
        (*context) -> calculateR_CPU(startLoc, startTime);
        (*context) -> calculateI_givenR_CPU(startLoc, startTime);
        ((*context) -> calculateP_SE_CPU(startLoc, startTime));
        return(0);
    }

    void FC_R_Star::sample(int verbose)
    {
        if (verbose){lssCout << "Sampling R_star\n";}
        (*currentSampler) -> drawSample();
    }

    long double FC_R_Star::getValue()
    {
        return(*(this -> value));
    }
    void FC_R_Star::setValue(long double val)
    {
        *(this -> value) = val;
    }

}
