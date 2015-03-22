#ifndef SPATIALSEIR_INITVALUE_CONTAINER
#define SPATIALSEIR_INITVALUE_CONTAINER
#include <Rcpp.h>
#include<ModelContext.hpp>
#include<modelComponent.hpp>

using namespace Rcpp;
using namespace SpatialSEIR;

RCPP_EXPOSED_CLASS(initialValueContainer)
class initialValueContainer : public modelComponent
{
    public:
        initialValueContainer();
        void setInitialValues(SEXP S0, SEXP E0, SEXP I0, SEXP R0,
                              SEXP S_star, SEXP E_star, SEXP I_star, SEXP R_star,
                              SEXP N);
        int getModelComponentType();
        int* compMatDim;
        int* S0;
        int* E0;
        int* I0;
        int* R0;
        int* S_star;
        int* E_star;
        int* I_star;
        int* R_star;
        int* N;
        ~initialValueContainer();
};

#endif
