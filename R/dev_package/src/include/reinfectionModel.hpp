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
        reinfectionModel(SEXP X, SEXP reinfectionMode);
        virtual void summary();
        int* xDim;
        double* X;
        double* betaPriorPrecision;
        ~reinfectionModel();
};

#endif
