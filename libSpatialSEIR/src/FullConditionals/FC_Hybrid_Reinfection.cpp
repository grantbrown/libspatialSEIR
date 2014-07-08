#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_FC_Hybrid_Reinfection.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    FC_Hybrid_Reinfection::FC_Hybrid_Reinfection(ModelContext * _context,
                         CompartmentalModelMatrix *_S_star, 
                         CompartmentalModelMatrix *_S, 
                         CompartmentalModelMatrix *_R, 
                         CompartmentalModelMatrix *_E_star,
                         CompartmentalModelMatrix *_R_star,
                         CompartmentFullConditional *_S_star_fc,
                         ParameterFullConditional *_beta_p_rs_fc,
                         InitData *_A0,
                         CovariateMatrix *_X, 
                         CovariateMatrix *_X_p_rs, 
                         double *_p_se,
                         double *_p_rs,
                         double *_beta_p_rs,
                         double _tausq,
                         double *_beta,
                         double *_rho,
                         double _steadyStateConstraintPrecision,
                         int _useOCL)
    {
       context = new ModelContext*;
       S_star = new CompartmentalModelMatrix*;
       S = new CompartmentalModelMatrix*;
       R = new CompartmentalModelMatrix*;
       E_star = new CompartmentalModelMatrix*;
       R_star = new CompartmentalModelMatrix*;
       parameterFullConditional = new ParameterFullConditional*;
       compartmentFullConditional = new CompartmentFullConditional*;
       A0 = new InitData*;
       X = new CovariateMatrix*;
       X_p_rs = new CovariateMatrix*;
       p_se = new double*;
       p_rs = new double*;
       beta_p_rs = new double*;
       tausq = new double;
       beta = new double*;
       rho = new double*;
       value = new long double;
       steadyStateConstraintPrecision = new double;
       useOCL = new int;
       samples = new int;
       accepted = new int;

       *context = _context;
       *S_star = _S_star;
       *S = _S;
       *R = _R;
       *E_star = _E_star;
       *R_star = _R_star;
       *parameterFullConditional = _beta_p_rs_fc; 
       *compartmentFullConditional = _S_star_fc;
       *A0 = _A0;
       *X = _X;
       *X_p_rs = _X_p_rs;
       *tausq = _tausq;
       *beta_p_rs = _beta_p_rs;
       *p_se = _p_se;
       *p_rs = _p_rs;
       *beta = _beta;
       *rho = _rho;
       *steadyStateConstraintPrecision = _steadyStateConstraintPrecision;
       *value = -1.0;
       *useOCL = _useOCL;
       *samples = 0;
       *accepted = 0;
    }    

    FC_Hybrid_Reinfection::~FC_Hybrid_Reinfection()
    {
        delete compartmentFullConditional;
        delete parameterFullConditional;
        delete S_star;
        delete S;
        delete R;
        delete E_star;
        delete R_star;
        delete A0;
        delete X;
        delete X_p_rs;
        delete p_se;
        delete p_rs;
        delete beta_p_rs;
        delete tausq;
        delete beta;
        delete rho;
        delete value;
        delete steadyStateConstraintPrecision;
        delete context;
        delete samples;
        delete accepted;
        delete useOCL;

    }

    int FC_Hybrid_Reinfection::evalCPU()
    {
        int nbeta = *((*X_p_rs) -> ncol_x);
        int nTpts = *((*R) -> nrow);
        int nLoc = *((*S)->ncol);
        int i,j,compIdx,Sstar_val,Estar_val,S_val,R_val;
        double p_se_val, p_rs_val;
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

        for (j = 0; j < nbeta; j++)
        {
            output -= ((*tausq)/2)*pow((*beta_p_rs)[j],2);
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

    int FC_Hybrid_Reinfection::evalOCL()
    {
        // Not Implemented
        return(evalCPU());
    }
    int FC_Hybrid_Reinfection::calculateRelevantCompartments()
    {
        (*context) -> calculateP_RS_CPU();
        (*context) -> calculateS_CPU();
        (*context) -> calculateR_givenS_CPU();
        return(0);
    }
    int FC_Hybrid_Reinfection::calculateRelevantCompartments_OCL()
    {
        (*context) -> calculateP_RS_CPU();
        (*context) -> calculateS_CPU();
        (*context) -> calculateR_givenS_CPU();
        return(0);
    }

    void FC_Hybrid_Reinfection::sample(int verbose)
    {
        if (verbose){std::cout << "Sampling Reinfection Parameters\n";}
        if (*useOCL){sampleOCL(); return;}
        sampleCPU();
    }

    int FC_Hybrid_Reinfection::sampleCPU()
    {
        this -> sampleHybrid_CPU(*context, 
                                 *beta_p_rs,
                                 *((*X_p_rs) -> ncol_x), 
                                 *S_star,
                                 *((*parameterFullConditional) -> sliceWidth),
                                 *((*compartmentFullConditional) -> sliceWidth));
        return 0;
    }
    int FC_Hybrid_Reinfection::sampleOCL()
    {
        // Not implemented
        this -> sampleCPU();
        return 0;
    }
    long double FC_Hybrid_Reinfection::getValue()
    {
        return(*(this -> value));
    }
    void FC_Hybrid_Reinfection::setValue(long double val)
    {
        *(this -> value) = val;
    }
}
