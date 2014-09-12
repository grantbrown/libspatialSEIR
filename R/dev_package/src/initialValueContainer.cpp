#include <Rcpp.h>
#include <initialValueContainer.hpp>
#include <DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

initialValueContainer::initialValueContainer()
{
    // Do nothing
}
void initialValueContainer::setInitialValues(SEXP S0_, SEXP E0_, SEXP I0_, SEXP R0_,
                                             SEXP S_star_, SEXP E_star_, SEXP I_star_,SEXP R_star_, SEXP N_)
{

    Rcpp::IntegerVector S0_vec(S0_);
    Rcpp::IntegerVector E0_vec(E0_);
    Rcpp::IntegerVector I0_vec(I0_);
    Rcpp::IntegerVector R0_vec(R0_);
    Rcpp::IntegerMatrix S_star_mat(S_star_);
    Rcpp::IntegerMatrix E_star_mat(E_star_);
    Rcpp::IntegerMatrix I_star_mat(I_star_);
    Rcpp::IntegerMatrix R_star_mat(R_star_);
    Rcpp::IntegerMatrix N_mat(N_);

    if (S0_vec.length() != E0_vec.length() || 
        E0_vec.length() != I0_vec.length() ||
        I0_vec.length() != R0_vec.length())
    {
        Rcpp::Rcout << "Init compartment lengths do not match.\n";
        throw(-1);
    }
    if (S_star_mat.nrow() != E_star_mat.nrow() ||
        E_star_mat.nrow() != I_star_mat.nrow() ||
        I_star_mat.nrow() != R_star_mat.nrow() || 
        R_star_mat.nrow() != N_mat.nrow()      ||
        S_star_mat.ncol() != E_star_mat.ncol() ||
        E_star_mat.ncol() != I_star_mat.ncol() ||
        I_star_mat.ncol() != R_star_mat.ncol() ||
        R_star_mat.ncol() != N_mat.ncol())
    {
        Rcpp::Rcout << "Compartment dimensions do not match\n";
        throw(-1);
    }

    int nTpt = (S_star_mat.nrow());
    int nLoc = (S_star_mat.ncol());
    compMatDim = new int[2];
    compMatDim[0] = nTpt;
    compMatDim[1] = nLoc;
    S0 = new int[nLoc];
    E0 = new int[nLoc];
    I0 = new int[nLoc];
    R0 = new int[nLoc];

    S_star = new int[nLoc*nTpt];
    E_star = new int[nLoc*nTpt];
    I_star = new int[nLoc*nTpt];
    R_star = new int[nLoc*nTpt];
    N = new int[nLoc*nTpt];

    memcpy(S0, S0_vec.begin(), nLoc*sizeof(int));
    memcpy(E0, E0_vec.begin(), nLoc*sizeof(int));
    memcpy(I0, I0_vec.begin(), nLoc*sizeof(int));
    memcpy(R0, R0_vec.begin(), nLoc*sizeof(int));
    memcpy(S_star, S_star_mat.begin(), nLoc*nTpt*sizeof(int));
    memcpy(E_star, E_star_mat.begin(), nLoc*nTpt*sizeof(int));
    memcpy(I_star, I_star_mat.begin(), nLoc*nTpt*sizeof(int));
    memcpy(R_star, R_star_mat.begin(), nLoc*nTpt*sizeof(int));
    memcpy(N, N_mat.begin(), nLoc*nTpt*sizeof(int));

}


initialValueContainer::~initialValueContainer()
{
    delete[] S0;
    delete[] E0;
    delete[] I0;
    delete[] R0;
    delete[] S_star;
    delete[] E_star;
    delete[] I_star;
    delete[] R_star;
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


