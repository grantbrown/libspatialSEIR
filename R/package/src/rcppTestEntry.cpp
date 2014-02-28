#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List rcppTestEntry()
{
    CharacterVector x = CharacterVector::create("Test1", "Test2");
    NumericVector y = NumericVector::create(0.0, 1.0);
    List z = List::create(x,y);
    return z;
}
