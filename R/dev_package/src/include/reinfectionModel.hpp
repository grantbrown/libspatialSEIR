#ifndef SPATIALSEIR_REINFECTION_MODEL
#define SPATIALSEIR_REINFECTION_MODEL
#include <Rcpp.h>
#include<ModelContext.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

RCPP_EXPOSED_CLASS(reinfectionModel)
class reinfectionModel
{
    public:
        reinfectionModel(SEXP reinfectionMode);
        virtual void buildReinfectionModel(SEXP _X, SEXP _paramInit, SEXP _prec);
        virtual void buildDummyReinfectionModel(int nTpt);
        virtual void summary();
        int* xDim;
        int* reinfectionMode;
        double* X;
        double* beta;
        double* betaPriorPrecision;
        ~reinfectionModel();
};

#endif
