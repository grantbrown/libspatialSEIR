#ifndef SPATIALSEIR_INITVALUE_CONTAINER
#define SPATIALSEIR_INITVALUE_CONTAINER
#include <Rcpp.h>
#include<ModelContext.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

RCPP_EXPOSED_CLASS(initialValueContainer)
class initialValueContainer
{
    public:
        initialValueContainer();
        void setInitialValues(SEXP S0, SEXP E0, SEXP I0, SEXP R0,
                              SEXP N);
        int* compMatDim;
        int* S0;
        int* E0;
        int* I0;
        int* R0;
        int* N;
        ~initialValueContainer();
};

#endif
