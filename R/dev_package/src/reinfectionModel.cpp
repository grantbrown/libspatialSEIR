#include <Rcpp.h>
#include <reinfectionModel.hpp>
#include <DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

reinfectionModel::reinfectionModel(SEXP _X, SEXP _prec)
{
    Rcpp::NumericMatrix inX(_X);
    Rcpp::NumericVector inPrecision(_prec);

    betaPriorPrecision = new double;
    *betaPriorPrecision = *(inPrecision.begin()); 
    xDim = new int[2];
    xDim[0] = inX.nrow();
    xDim[1] = inX.ncol();
    X = new double[xDim[0]*xDim[1]];
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
    delete[] X;
    delete betaPriorPrecision;
}

RCPP_MODULE(mod_reinfectionModel)
{
    using namespace Rcpp;
    class_<reinfectionModel>( "reinfectionModel" )
    .constructor<SEXP,SEXP>()
    .method("summary", &reinfectionModel::summary);
}


