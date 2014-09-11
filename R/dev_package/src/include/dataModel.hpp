#ifndef SPATIALSEIR_DATA_MODEL
#define SPATIALSEIR_DATA_MODEL
#include <Rcpp.h>
#include<ModelContext.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

RCPP_EXPOSED_CLASS(dataModel)
class dataModel
{
    public:
        dataModel(SEXP Y, SEXP type);
        virtual void summary();
        int* nLoc;
        int* nTpt;
        int* dataModelType;
        int* Y;

        ~dataModel();
};

#endif
