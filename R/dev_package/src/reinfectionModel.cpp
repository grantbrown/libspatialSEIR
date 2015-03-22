#include <Rcpp.h>
#include <reinfectionModel.hpp>
#include <DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;


reinfectionModel::reinfectionModel(SEXP reinfectMode)
{
     Rcpp::IntegerVector modeVec(reinfectMode);
     reinfectionMode = new int;
     *reinfectionMode = modeVec[0];
     betaPriorPrecision = new double; *betaPriorPrecision = -1; 
     xDim = new int[2]; memset(xDim, 1, 2*sizeof(int));
}

int reinfectionModel::getModelComponentType()
{
    return(LSS_REINFECTION_MODEL_TYPE);
}

void reinfectionModel::buildReinfectionModel(SEXP _X, SEXP _paramInit, SEXP _priorMean, SEXP _prec)
{
    Rcpp::NumericMatrix inX(_X);
    Rcpp::NumericVector inPrecision(_prec);
    Rcpp::NumericVector initValues(_paramInit);
    Rcpp::NumericVector priorMeans(_priorMean);
    xDim[0] = inX.nrow();
    xDim[1] = inX.ncol();
    beta = new double[xDim[1]];
    betaPriorMean = new double[xDim[1]];
    betaPriorPrecision = new double[xDim[1]];
    X = new double[xDim[0]*xDim[1]];
    if (initValues.length() != xDim[1] || 
        priorMeans.length() != xDim[1] || 
        inPrecision.length() != xDim[1])
    {
        ::Rf_error("Number of parameters, prior means, or precisions does not equal the number of supplied covariates.\n");
    }
    int i;
    for (i = 0; i < xDim[1]; i++)
    {
        beta[i] = initValues[i];
        betaPriorPrecision[i] = inPrecision[i];
        betaPriorMean[i] = priorMeans[i];
    }
    for (i = 0; i < xDim[0]*xDim[1]; i++)
    {
       X[i] = inX[i]; 
    }
}

void reinfectionModel::buildDummyReinfectionModel(int nTpt)
{
    xDim[0] = nTpt;
    xDim[1] = 1;
    X = new double[nTpt];
    int i;
    for (i = 0; i < nTpt; i++)
    {
        X[i] = 1.0;
    }
    beta = new double[1];
    betaPriorMean = new double[1];
    betaPriorPrecision = new double[1];
    *beta = -INFINITY;
    *betaPriorPrecision = -INFINITY;
    *betaPriorMean = -INFINITY;
    *betaPriorPrecision=100;
}

void reinfectionModel::summary()
{
    Rcpp::Rcout << "Reinfection covariate dimensions: " << xDim[0] << ", " << xDim[1] << "\n";
    Rcpp::Rcout << "Linear model parameter prior precision: ";
    int i;
    for (i = 0; i < xDim[1] - 1; i++)
    {
        Rcpp::Rcout << betaPriorPrecision[i] << ", "; 
    }
    Rcpp::Rcout << "Linear model parameter prior means: ";
    for (i = 0; i < xDim[1] - 1; i++)
    {
        Rcpp::Rcout << betaPriorMean[i] << ", "; 
    }
    Rcpp::Rcout << betaPriorMean[xDim[1]-1] << "\n";
}
reinfectionModel::~reinfectionModel()
{
    delete[] xDim;
    delete[] beta;
    delete[] X;
    delete[] betaPriorPrecision;
    delete[] betaPriorMean;
}

RCPP_MODULE(mod_reinfectionModel)
{
    using namespace Rcpp;
    class_<reinfectionModel>( "reinfectionModel" )
    .constructor<SEXP>()
    .method("buildReinfectionModel", &reinfectionModel::buildReinfectionModel)
    .method("summary", &reinfectionModel::summary);
}


