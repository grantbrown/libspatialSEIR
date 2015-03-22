#ifndef SPATIALSEIR_DATA_MODEL
#define SPATIALSEIR_DATA_MODEL
#include <Rcpp.h>
#include<ModelContext.hpp>
#include<modelComponent.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

RCPP_EXPOSED_CLASS(dataModel)
class dataModel : public modelComponent
{
    public:
        dataModel(SEXP Y, SEXP type, SEXP compartment);
        virtual void summary();
        virtual void setOverdispersionParameters(SEXP priorAlpha, SEXP priorBeta, SEXP initialValue);
        Rcpp::IntegerVector* compartmentDimensions;
        int getModelComponentType();
        double* priorParameters;
        double* initialParameterValues;
        int* nLoc;
        int* nTpt;
        int* dataModelType;
        int* dataModelCompartment;
        int* Y;
        int* setMode;

        ~dataModel();
};

#endif
