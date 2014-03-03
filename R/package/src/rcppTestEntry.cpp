#include <Rcpp.h>
#include <CompartmentalModelMatrix.hpp>
#include <CovariateMatrix.hpp>
#include <OCLProvider.hpp>

using namespace Rcpp;
using namespace SpatialSEIR;
// [[Rcpp::export]]
List rcppTestEntry()
{
    CovariateMatrix *CovMat = new CovariateMatrix();
    CompartmentalModelMatrix *CompMat = new CompartmentalModelMatrix();
    OCLProvider *OCLObj = new OCLProvider();

    NumericVector y = NumericVector::create(0.0, 1.0);
    List z = List::create(y);
    return z;
}
