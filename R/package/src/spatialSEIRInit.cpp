#include <Rcpp.h>
#include <FullConditional.hpp>

using namespace Rcpp;
using namespace SpatialSEIR;
// [[Rcpp::export]]
List spatialSEIRInit(SEXP compMatDim,
                     SEXP xDim,
                     SEXP S0_,
                     SEXP E0_,
                     SEXP I0_,
                     SEXP R0_,
                     SEXP Sstar0,
                     SEXP Estar0,
                     SEXP Istar0,
                     SEXP Rstar0,
                     SEXP Istar, 
                     SEXP X_)
{

    //Deal with the data conversion from R to c++
    Rcpp::IntegerVector compartmentDimensions(compMatDim);
    Rcpp::IntegerVector covariateDimensions(xDim);
    Rcpp::IntegerVector S0(S0_);
    Rcpp::IntegerVector E0(E0_);
    Rcpp::IntegerVector I0(I0_);
    Rcpp::IntegerVector R0(R0_);

    Rcpp::IntegerVector S_star0(Sstar0);
    Rcpp::IntegerVector E_star0(Estar0);
    Rcpp::IntegerVector I_star0(Istar0);
    Rcpp::IntegerVector R_star0(Rstar0);

    Rcpp::IntegerVector I_star(Istar);
    Rcpp::NumericVector X(X_);


    CovariateMatrix *CovMat = new CovariateMatrix();


    CovMat -> genFromDataStream(X.begin(), 
                                &covariateDimensions[0], 
                                &covariateDimensions[1]);


    CompartmentalModelMatrix *CompMat = new CompartmentalModelMatrix();
    CompMat -> genFromDataStream(I_star.begin(), 
                                 &compartmentDimensions[0],
                                 &compartmentDimensions[1]);

    Rcpp::IntegerVector y = Rcpp::IntegerVector::create(0);
    List z = List::create(y);
    return z;
}
