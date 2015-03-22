#ifndef SPATIALSEIR_EXPOSURE_MODEL
#define SPATIALSEIR_EXPOSURE_MODEL
#include <Rcpp.h>
#include<ModelContext.hpp>
#include<modelComponent.hpp>

using namespace Rcpp;
using namespace SpatialSEIR;

RCPP_EXPOSED_CLASS(exposureModel)
class exposureModel : public modelComponent
{
    public:
        exposureModel(SEXP X, SEXP ntpt, SEXP nloc, SEXP initBeta, SEXP priorMean, SEXP precision);
        virtual void summary();
        int getModelComponentType();
        virtual Rcpp::NumericVector getOffset();
        virtual void setOffset(Rcpp::NumericVector offs);
        double* offset;
        int* xDim;
        int* nTpt;
        int* nLoc;
        double* beta;
        double* X;
        double* betaPriorPrecision;
        double* betaPriorMean;
        ~exposureModel();
};

#endif
