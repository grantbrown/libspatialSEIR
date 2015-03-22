#include <Rcpp.h>
#include <distanceModel.hpp>
#include <DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

distanceModel::distanceModel()
{
    numLocations = new int; *numLocations = -1;
    scaledDistArgs = new scaledDistanceArgs();    
    scaledDistArgs -> priorAlpha_rho = 1.0;
    scaledDistArgs -> priorBeta_rho = 1.0;
}

int distanceModel::getModelComponentType()
{
    return(LSS_DISTANCE_MODEL_TYPE);
}

void distanceModel::setPriorParameters(SEXP arg1, SEXP arg2)
{
    Rcpp::NumericVector a(arg1);
    Rcpp::NumericVector b(arg2);
    scaledDistArgs -> priorAlpha_rho = a[0];
    scaledDistArgs -> priorBeta_rho = b[0];
    if (a[0] < 0 || b[0] < 0)
    {
        Rcpp::Rcout << "Invalid prior values: " << a[0] << ", " << b[0] << ", " << "setting to uniform prior\n";
        scaledDistArgs -> priorAlpha_rho = 1.0;
        scaledDistArgs -> priorBeta_rho = 1.0;
    }
}

void distanceModel::addDistanceMatrix(NumericMatrix distMat)
{
    if (distMat.nrow() != distMat.ncol())
    {
        ::Rf_error("Distance matrix must be square.\n");
    }
    else if (*numLocations != -1 && distMat.nrow() != (*numLocations))
    {
        ::Rf_error("Dimension does not match previously added distance matrix\n");
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
    .method("setPriorParameters", &distanceModel::setPriorParameters)
    .property("numMatrices", &distanceModel::getNumDistanceMatrices, "Number of distict distance matrices.");
}


