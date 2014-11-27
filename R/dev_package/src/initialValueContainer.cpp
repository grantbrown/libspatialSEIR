#include <Rcpp.h>
#include <initialValueContainer.hpp>
#include <DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

initialValueContainer::initialValueContainer()
{
    // Do nothing
}
void initialValueContainer::setInitialValues(SEXP S0_, SEXP E0_, SEXP I0_, SEXP R0_, SEXP N_)
{

    Rcpp::IntegerVector S0_vec(S0_);
    Rcpp::IntegerVector E0_vec(E0_);
    Rcpp::IntegerVector I0_vec(I0_);
    Rcpp::IntegerVector R0_vec(R0_);

    Rcpp::IntegerMatrix N_mat(N_);

    if (S0_vec.length() != E0_vec.length() || 
        E0_vec.length() != I0_vec.length() ||
        I0_vec.length() != R0_vec.length())
    {
        Rcpp::Rcout << "Init compartment lengths do not match.\n";
        throw(-1);
    }

    compMatDim = new int[2];
    int nTpt = N_mat.nrow();
    int nLoc = N_mat.ncol();
    compMatDim[0] = nTpt;
    compMatDim[1] = nLoc;

    S0 = new int[nLoc];
    E0 = new int[nLoc];
    I0 = new int[nLoc];
    R0 = new int[nLoc];

    N = new int[nLoc*nTpt];

    int i;
    for (i = 0; i < nLoc; i++)
    {
        S0[i] = S0_vec[i];
        E0[i] = E0_vec[i];
        I0[i] = I0_vec[i];
        R0[i] = R0_vec[i];
    }
    for (i = 0; i < nLoc*nTpt; i++)
    {
        N[i] = N_mat[i];
    }
}


initialValueContainer::~initialValueContainer()
{
    delete[] S0;
    delete[] E0;
    delete[] I0;
    delete[] R0;
    delete[] compMatDim;
    delete[] N;
}

RCPP_MODULE(mod_initialValueContainer)
{
    using namespace Rcpp;
    class_<initialValueContainer>( "initialValueContainer" )
    .constructor()
    .method("setInitialValues", &initialValueContainer::setInitialValues);
}


