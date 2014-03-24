#include <Rcpp.h>
#include <cmath>
#include <ModelContext.hpp>

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
                     SEXP Sstar, 
                     SEXP Estar, 
                     SEXP Istar, 
                     SEXP Rstar, 
                     SEXP X_,
                     SEXP Z_,
                     SEXP DistMat_,
                     SEXP rho_,
                     SEXP beta_,
                     SEXP p_ei_,
                     SEXP p_ir_,
                     SEXP p_rs_,
                     SEXP N_)
{
    Rcpp::Rcout << "Wrapping input data in Rcpp vectors.\n";
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

    Rcpp::IntegerVector S_star(Sstar);
    Rcpp::IntegerVector E_star(Estar);
    Rcpp::IntegerVector I_star(Istar);
    Rcpp::IntegerVector R_star(Rstar);

    Rcpp::NumericVector X(X_);
    Rcpp::NumericVector Z(Z_);
    Rcpp::NumericVector DistMat(DistMat_);

    Rcpp::NumericVector rho(rho_);
    Rcpp::NumericVector beta(beta_);
    Rcpp::NumericVector p_ei(p_ei_);
    Rcpp::NumericVector p_ir(p_ir_);
    Rcpp::NumericVector p_rs(p_rs_);
    Rcpp::IntegerVector N(N_);

    Rcpp::Rcout << "Creating Model Context\n";
    // Create the empty ModelContext object  
    ModelContext* context = new ModelContext();

    Rcpp::Rcout << "Loading covariate information into model context object\n";
    // Create the covariate matrix object. 
    context -> X -> genFromDataStream(X.begin(), 
                                Z.begin(),
                                &covariateDimensions_x[0], 
                                &covariateDimensions_x[1],
                                &covariateDimensions_z[0], 
                                &covariateDimensions_z[1]);

    Rcpp::Rcout << "Creating empty S,E,I,R matrices.\n";
    // Populate the CompartmentalModelMatrix objects.  
    context -> S -> createEmptyCompartment(&compartmentDimensions[0],
                                           &compartmentDimensions[1]); 
    context -> E -> createEmptyCompartment(&compartmentDimensions[0],
                                           &compartmentDimensions[1]); 
    context -> I -> createEmptyCompartment(&compartmentDimensions[0],
                                           &compartmentDimensions[1]); 
    context -> R -> createEmptyCompartment(&compartmentDimensions[0],
                                           &compartmentDimensions[1]); 

    Rcpp::Rcout << "Filling initial S_star, E_star, I_star_ R_star matrices.\n";
    context -> S_star -> genFromDataStream(S_star.begin(), 
                                           &compartmentDimensions[0],
                                           &compartmentDimensions[1]);
    context -> E_star -> genFromDataStream(E_star.begin(), 
                                           &compartmentDimensions[0],
                                           &compartmentDimensions[1]);
    context -> I_star -> genFromDataStream(I_star.begin(), 
                                           &compartmentDimensions[0],
                                           &compartmentDimensions[1]);
    context -> R_star -> genFromDataStream(R_star.begin(), 
                                           &compartmentDimensions[0],
                                           &compartmentDimensions[1]);

    Rcpp::Rcout << "Rcpp Provided Num Locations: " << compartmentDimensions[0] << "\n";
    Rcpp::Rcout << "Rcpp Provided Num Times: " << compartmentDimensions[1] << "\n";

    Rcpp::Rcout << "Creating raw and scaled distance matrices.\n";
    context -> rawDistMat -> genFromDataStream(DistMat.begin(), &compartmentDimensions[0]);
    context -> scaledDistMat -> genFromDataStream(DistMat.begin(), &compartmentDimensions[0]);
    context -> scaledDistMat -> scaledInvFunc_CPU(60*60*2, context -> rawDistMat -> data);

    // Populate the Time 0 initialization data
    context -> A0 ->  populate(S0.begin(),E0.begin(),I0.begin(),R0.begin(),
                               S_star0.begin(),E_star0.begin(),I_star0.begin(),
                               R_star0.begin(),&compartmentDimensions[0]);
    Rcpp::Rcout << "Populating Model Context\n";
    context -> populate(rho.begin(), beta.begin(), p_ei.begin(), 
                        p_ir.begin(),p_rs.begin(),N.begin());
    Rcpp::Rcout << "Calculating Eta\n";
    context -> X -> calculate_eta_CPU(context -> eta, context -> beta);

    // Test calculation functions. 
    
    Rcpp::Rcout << "Calculating S\n";
    context -> calculateS_CPU();

    Rcpp::Rcout << "Calculating E\n";
    context -> calculateE_CPU();
    
    Rcpp::Rcout << "Calculating I\n";
    context -> calculateI_CPU();

    Rcpp::Rcout << "Calculating R\n";
    context -> calculateR_CPU();

    Rcpp::Rcout << "Calculating P_SE from Starting Beta, I\n";
    context -> calculateP_SE_CPU();

    Rcpp::Rcout << "Testing S_star FC: \n";
    int tmp = context -> S_star_fc -> evalCPU();
    Rcpp::Rcout << "returned: " << tmp << ", value: " << *(context -> S_star_fc -> value) << "\n";
    Rcpp::Rcout << "Testing E_star FC: \n";
    tmp = context -> E_star_fc -> evalCPU();
    Rcpp::Rcout << "returned: " << tmp << ", value: " << *(context -> E_star_fc -> value) << "\n";
    Rcpp::Rcout << "Testing R_star FC: \n";
    tmp = context -> R_star_fc -> evalCPU();
    Rcpp::Rcout << "returned: " << tmp << ", value: " << *(context -> R_star_fc -> value) << "\n";
    Rcpp::Rcout << "Testing beta FC: \n";
    tmp = context -> beta_fc -> evalCPU();
    Rcpp::Rcout << "returned: " << tmp << ", value: " << *(context -> beta_fc -> value) << "\n";
    Rcpp::Rcout << "Testing p_rs FC: \n";
    tmp = context -> p_rs_fc -> evalCPU();
    Rcpp::Rcout << "returned: " << tmp << ", value: " << *(context -> p_rs_fc -> value) << "\n";
    Rcpp::Rcout << "Testing rho FC: \n";
    tmp = context -> rho_fc -> evalCPU();
    Rcpp::Rcout << "returned: " << tmp << ", value: " << *(context -> rho_fc -> value) << "\n";
    Rcpp::Rcout << "Testing p_ei FC: \n";
    tmp = context -> p_ei_fc -> evalCPU();
    Rcpp::Rcout << "returned: " << tmp << ", value: " << *(context -> p_ei_fc -> value) << "\n";

    Rcpp::XPtr<ModelContext*> ptr(&context, true);

    // Clean up
    Rcpp::Rcout << "Finished.\n";
    return ptr;
}
