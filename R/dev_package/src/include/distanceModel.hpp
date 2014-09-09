#ifndef SPATIALSEIR_DISTANCE_MODEL
#define SPATIALSEIR_DISTANCE_MODEL
#include <Rcpp.h>
#include<DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

class distanceModel
{
    public:
        distanceModel();
        virtual void addDistanceMatrix(NumericMatrix distMat);
        virtual void summary();

        int* numLocations;
        scaledDistanceArgs* scaledDistArgs;

        ~distanceModel();
};

#endif
