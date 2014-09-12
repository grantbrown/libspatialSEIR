#include <Rcpp.h>
#include <reinfectionModel.hpp>
#include <DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;


reinfectionModel::reinfectionModel(SEXP reinfectMode)
{
     Rcpp::Integervector modeVec(reinfectMode);
     reinfectionMode = new int;
     *reinfectionMode = modeVec[0];
     betaPriorPrecision = new double; *betaPriorPrecision = -1; 
     xDim = new int[2]; memset(xDim, -1, 2*sizeof(int));
}

reinfectionModel::buildReinfectionModel(SEXP _X, SEXP _paramInit, SEXP _prec)
{
    Rcpp::NumericMatrix inX(_X);
    Rcpp::NumericVector inPrecision(_prec);
    Rcpp::NumericVector initValues(_paramInit);
    *betaPriorPrecision = *(inPrecision.begin()); 
    xDim[0] = inX.nrow();
    xDim[1] = inX.ncol();
    beta = new double[xDim[1]];
    X = new double[xDim[0]*xDim[1]];
    if (initValues.length() != xDim[1])
    {
        Rcpp::Rcout << "Number of parameters does not equal the number of supplied covariates.\n";
        throw(-1);
    }
    memcpy(beta, initValues.begin(), xDim[1]*sizeof(double));
    memcpy(X, inX.begin(), xDim[0]*xDim[1]*sizeof(double));
}

void reinfectionModel::summary()
{
    Rcpp::Rcout << "Reinfection covariate dimensions: " << xDim[0] << ", " << xDim[1] << "\n";
    Rcpp::Rcout << "Linear model parameter prior precision: " << *betaPriorPrecision << "\n";

}
reinfectionModel::~reinfectionModel()
{
    delete[] xDim;
    delete[] beta;
    delete[] X;
    delete betaPriorPrecision;
}

RCPP_MODULE(mod_reinfectionModel)
{
    using namespace Rcpp;
    class_<reinfectionModel>( "reinfectionModel" )
    .constructor()
    .method("buildReinfectionModel", &reinfectionModel::buildReinfectionModel)
    .method("summary", &reinfectionModel::summary);
}


