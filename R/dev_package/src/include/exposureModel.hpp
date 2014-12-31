#ifndef SPATIALSEIR_EXPOSURE_MODEL
#define SPATIALSEIR_EXPOSURE_MODEL
#include <Rcpp.h>
#include<ModelContext.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

RCPP_EXPOSED_CLASS(exposureModel)
class exposureModel
{
    public:
        exposureModel(SEXP X, SEXP Z, SEXP initBeta, SEXP priorMean, SEXP precision, SEXP hasZ);
        virtual void summary();
        virtual Rcpp::NumericVector getOffset();
        virtual void setOffset(Rcpp::NumericVector offs);
        double* offset;
        int* xDim;
        int* zDim;
        double* beta;
        double* X;
        double* Z;
        double* betaPriorPrecision;
        double* betaPriorMean;
        ~exposureModel();
};

#endif
