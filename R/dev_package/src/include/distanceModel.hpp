#ifndef SPATIALSEIR_DISTANCE_MODEL
#define SPATIALSEIR_DISTANCE_MODEL
#include <Rcpp.h>
#include<DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

RCPP_EXPOSED_CLASS(distanceModel)
class distanceModel
{
    public:
        distanceModel();
        virtual void addDistanceMatrix(NumericMatrix distMat);
        virtual void setPriorParameters(SEXP a, SEXP b);
        virtual void summary();
        virtual int getNumDistanceMatrices();

        int* numLocations;
        double* priorAlpha;
        double* priorBeta;
        scaledDistanceArgs* scaledDistArgs;

        ~distanceModel();
};

#endif
