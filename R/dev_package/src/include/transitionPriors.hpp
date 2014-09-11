#ifndef SPATIALSEIR_TRANSITION_PRIORS
#define SPATIALSEIR_TRANSITION_PRIORS
#include <Rcpp.h>
#include<ModelContext.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

RCPP_EXPOSED_CLASS(transitionPriors)
class transitionPriors
{
    public:
        transitionPriors();
        void setPriorsFromProbabilities(SEXP p_ei, SEXP p_ir, SEXP p_ei_ess, SEXP p_ir_ess);
        void setUniformPriors();
        void setPriorsManually(SEXP priorAlpha_gammaEI, SEXP priorBeta_gammaEI,
                               SEXP priorAlpha_gammaIR, SEXP priorBeta_gammaIR);
        void summary();

        double* gamma_ei;
        double* gamma_ir;
        double* gamma_ei_params;
        double* gamma_ir_params;
        ~transitionPriors();
};

#endif
