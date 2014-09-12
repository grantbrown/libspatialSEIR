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
        exposureModel(SEXP X, SEXP Z, SEXP initBeta, SEXP precision);
        virtual void summary();
        virtual Rcpp::IntegerVector getOffset();
        virtual void setOffset(Rcpp::IntegerVector offs);
        double* offset;
        int* xDim;
        int* zDim;
        double* beta;
        double* X;
        double* Z;
        double* betaPriorPrecision;
        ~exposureModel();
};

#endif
