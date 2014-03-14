#include <Rcpp.h>
#include <cmath>
#include <FullConditional.hpp>

using namespace Rcpp;
using namespace SpatialSEIR;
// [[Rcpp::export]]
SEXP spatialSEIRInit(SEXP compMatDim,
                     SEXP xDim,
                     SEXP zDim,
                     SEXP S0_,
                     SEXP E0_,
                     SEXP I0_,
                     SEXP R0_,
                     SEXP Sstar0,
                     SEXP Estar0,
                     SEXP Istar0,
                     SEXP Rstar0,
                     SEXP Istar, 
                     SEXP X_,
                     SEXP Z_)
{

    //Deal with the data conversion from R to c++
    Rcpp::IntegerVector compartmentDimensions(compMatDim);
    Rcpp::IntegerVector covariateDimensions_x(xDim);
    Rcpp::IntegerVector covariateDimensions_z(zDim);
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
    Rcpp::NumericVector Z(Z_);

    CovariateMatrix *CovMat = new CovariateMatrix();


    CovMat -> genFromDataStream(X.begin(), 
                                Z.begin(),
                                &covariateDimensions_x[0], 
                                &covariateDimensions_x[1],
                                &covariateDimensions_z[0], 
                                &covariateDimensions_z[1]);
    // Test CPU eta calculation. 

    int i;
    double* beta = new double[covariateDimensions_x[1]];
    double* gamma = new double[covariateDimensions_z[1]];
    double* eta = new double[covariateDimensions_x[1] + covariateDimensions_z[1]];


    for (i = 0; i < covariateDimensions_x[1]; i++)
    {
        beta[i] = 0.0;
    }


    for (i = 0; i < covariateDimensions_z[1]; i++)
    {
        gamma[i] = 0.0;
    }



    CovMat -> calculate_eta_CPU(eta, beta, gamma);
    Rcpp::Rcout << "Check eta\n";

    for (i = 0; i < (covariateDimensions_x[1] + covariateDimensions_z[1]); i++)
    {
        if (std::abs(eta[i]) > 0.000000001)
        {
            Rcpp::Rcout << "Eta " << i << " is not zero\n";
        }
    }



    CompartmentalModelMatrix *CompMat = new CompartmentalModelMatrix();
    CompMat -> genFromDataStream(I_star.begin(), 
                                 &compartmentDimensions[0],
                                 &compartmentDimensions[1]);




    Rcpp::IntegerVector y = Rcpp::IntegerVector::create(0);
    //List z = List::create(y);
    Rcpp::XPtr<CompartmentalModelMatrix*> ptr(&CompMat, true);


    // Clean up

    delete[] beta;
    delete[] gamma;
    delete[] eta;

    Rcpp::Rcout << "Finished.\n";
    return ptr;
}
