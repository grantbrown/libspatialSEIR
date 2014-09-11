#include <Rcpp.h>
#include <Rmath.h>
#include <transitionPriors.hpp>
#include <DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

transitionPriors::transitionPriors()
{
    gamma_ei = new double;
    gamma_ir = new double;
    gamma_ei_params = new double[2];
    gamma_ir_params = new double[2];
}

void transitionPriors::setUniformPriors()
{
    gamma_ei_params[0] = 1.0; 
    gamma_ei_params[1] = 1.0;
    gamma_ir_params[0] = 1.0; 
    gamma_ir_params[1] = 1.0;

    *gamma_ei = R::rgamma(gamma_ei_params[0], gamma_ei_params[1]);
    *gamma_ir = R::rgamma(gamma_ir_params[0], gamma_ir_params[1]);

}

void transitionPriors::setPriorsFromProbabilities(SEXP p_ei, SEXP p_ir, 
                                                  SEXP p_ei_ess, SEXP p_ir_ess)
{
    double pEI, pIR;
    int pEIess, pIRess;
    Rcpp::NumericVector p_ei_vec(p_ei);
    Rcpp::NumericVector p_ir_vec(p_ir);
    Rcpp::IntegerVector p_ei_ess_vec(p_ei_ess);
    Rcpp::IntegerVector p_ir_ess_vec(p_ir_ess);

    pEI = p_ei_vec[0]; pIR = p_ir_vec[0];
    pEIess = p_ei_ess_vec[0]; pIRess = p_ir_ess_vec[0];


    *gamma_ei = -std::log(1-pEI);
    *gamma_ir = -std::log(1-pIR);

    gamma_ei_params[0] = pEIess;
    gamma_ei_params[1] = pEIess/(*gamma_ei);

    gamma_ir_params[0] = pIRess;
    gamma_ir_params[1] = pIRess/(*gamma_ir);
     
}

void transitionPriors::setPriorsManually(SEXP priorAlpha_gammaEI, SEXP priorBeta_gammaEI,
                       SEXP priorAlpha_gammaIR, SEXP priorBeta_gammaIR)
{
    Rcpp::NumericVector pA_gammaEI(priorAlpha_gammaEI);
    Rcpp::NumericVector pB_gammaEI(priorBeta_gammaEI);
    Rcpp::NumericVector pA_gammaIR(priorAlpha_gammaIR);
    Rcpp::NumericVector pB_gammaIR(priorBeta_gammaIR);

    gamma_ei_params[0] = pA_gammaEI[0]; 
    gamma_ei_params[1] = pB_gammaEI[0]; 

    gamma_ir_params[0] = pA_gammaIR[0]; 
    gamma_ir_params[1] = pB_gammaIR[0]; 

    *gamma_ei = R::rgamma(gamma_ei_params[0], gamma_ei_params[1]);
    *gamma_ir = R::rgamma(gamma_ir_params[0], gamma_ir_params[1]);
}

void transitionPriors::summary()
{
    Rcpp::Rcout << "gamma_ei: " << *gamma_ei << "\n";
    Rcpp::Rcout << "gamma_ir: " << *gamma_ir << "\n";
    Rcpp::Rcout << "gamma_ei parameters: " << gamma_ei_params[0] << ", " << gamma_ei_params[1] << "\n";
    Rcpp::Rcout << "gamma_ir parameters: " << gamma_ir_params[0] << ", " << gamma_ir_params[1] << "\n";


}
transitionPriors::~transitionPriors()
{
    delete gamma_ei;
    delete gamma_ir;
    delete[] gamma_ei_params;
    delete[] gamma_ir_params;
}

RCPP_MODULE(mod_transitionPriors)
{
    using namespace Rcpp;
    class_<transitionPriors>( "transitionPriors" )
    .constructor()
    .method("setUniformPriors", &transitionPriors::setUniformPriors)
    .method("setPriorsFromProbabilities", &transitionPriors::setPriorsFromProbabilities)
    .method("setPriorsManually", &transitionPriors::setPriorsManually)
    .method("summary", &transitionPriors::summary);
}


