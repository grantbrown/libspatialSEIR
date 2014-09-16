#include <Rcpp.h>
#include <distanceModel.hpp>
#include <DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

distanceModel::distanceModel()
{
    numLocations = new int; *numLocations = -1;
    scaledDistArgs = new scaledDistanceArgs();    
}

void distanceModel::addDistanceMatrix(NumericMatrix distMat)
{
    if (distMat.nrow() != distMat.ncol())
    {
        Rcpp::Rcout << "Distance matrix must be square.\n";
        throw(-1);
    }
    else if (*numLocations != -1 && distMat.nrow() != (*numLocations))
    {
        Rcpp::Rcout << "Dimension does not match previously added distance matrix.\n";
        throw(-1);
    }
    double* distAlloc = new double[distMat.nrow()*distMat.ncol()];
    int i;
    for (i = 0; i < distMat.nrow()*distMat.ncol(); i++)
    {
        distAlloc[i] = distMat[i];
    }
    scaledDistArgs -> dim = numLocations;
    (scaledDistArgs -> inData).push_back(distAlloc);
    *numLocations = distMat.nrow();
}
void distanceModel::summary()
{
    Rcpp::Rcout << "Number of locations: " << *numLocations << "\n";
    Rcpp::Rcout << "Number of distance structures: " << ((scaledDistArgs -> inData).size()) << "\n";
}
distanceModel::~distanceModel()
{
    unsigned int k;
    unsigned int mats = (scaledDistArgs -> inData.size());
    for (k = 0; k < mats; k++)
    {
        delete[] ((scaledDistArgs -> inData)[k]);
    }
    delete scaledDistArgs;
    delete numLocations;
}

int distanceModel::getNumDistanceMatrices()
{
    return((scaledDistArgs -> inData).size()); 
}

RCPP_MODULE(mod_distanceModel)
{
    using namespace Rcpp;
    class_<distanceModel>( "distanceModel" )
    .constructor()
    .method("addDistanceMatrix", &distanceModel::addDistanceMatrix)
    .method("summary", &distanceModel::summary)
    .property("numMatrices", &distanceModel::getNumDistanceMatrices, "Number of distict distance matrices.");
}


