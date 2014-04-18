#include <Rcpp.h>
#include <cmath>
#include <ModelContext.hpp>
#include <FullConditional.hpp>
#include <CovariateMatrix.hpp>
#include <CompartmentalModelMatrix.hpp>
#include <DistanceMatrix.hpp>
#include <IOProvider.hpp>
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
                     SEXP gamma_,
                     SEXP priorAlpha_gamma_,
                     SEXP priorBeta_gamma_,
                     SEXP beta_,
                     SEXP p_ei_,
                     SEXP p_ir_,
                     SEXP p_rs_,
                     SEXP N_,
                     SEXP outFile,
                     SEXP logVarList,
                     SEXP iterationStride,
                     SEXP verboseFlag,
                     SEXP debugFlag)
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

    Rcpp::NumericVector gamma(gamma_);
    Rcpp::NumericVector priorAlpha_gamma(priorAlpha_gamma_);
    Rcpp::NumericVector priorBeta_gamma(priorBeta_gamma_);

    Rcpp::NumericVector beta(beta_);
    Rcpp::NumericVector p_ei(p_ei_);
    Rcpp::NumericVector p_ir(p_ir_);
    Rcpp::NumericVector p_rs(p_rs_);
    Rcpp::IntegerVector N(N_);

    std::string* chainOutputFile = new std::string(); 
    *chainOutputFile = Rcpp::as<std::string>(outFile);
    Rcpp::IntegerVector chainOutputControl(logVarList); 
    // logVarList: (Beta, rho,gamma,p_se,p_ei,p_ir,p_rs,S*, E*, I*, R*)
    // Nonzero if respective variables are to be output
    Rcpp::IntegerVector chainStride(iterationStride);
    
    Rcpp::IntegerVector verbose(verboseFlag);
    Rcpp::IntegerVector debug(debugFlag);

    // Sanity check the input data. 
    int compartmentSize = (compartmentDimensions[0]*compartmentDimensions[1]);
    if (S_star0.size() != compartmentDimensions[0])
    {
        Rcpp::Rcout << "Invalid S_star0 Compartment Size!\n";
        throw(-1);
    }
    if (E_star0.size() != compartmentDimensions[0])
    {
        Rcpp::Rcout << "Invalid E_star0 Compartment Size!\n";
        throw(-1);
    }
    if (I_star0.size() != compartmentDimensions[0])
    {
        Rcpp::Rcout << "Invalid I_star0 Compartment Size!\n";
        throw(-1);
    }
    if (R_star0.size() != compartmentDimensions[0])
    {
        Rcpp::Rcout << "Invalid R_star0 Compartment Size!\n";
        throw(-1);
    }


    if (S_star.size() != compartmentSize)
    {
        Rcpp::Rcout << "Invalid S_star Compartment Size!\n";
        throw(-1);
    }
    if (E_star.size() != compartmentSize)
    {
        Rcpp::Rcout << "Invalid E_star Compartment Size!\n";
        throw(-1);
    }
    if (I_star.size() != compartmentSize)
    {
        Rcpp::Rcout << "Invalid I_star Compartment Size!\n";
        throw(-1);
    }
    if (R_star.size() != compartmentSize)
    {
        Rcpp::Rcout << "Invalid R_star Compartment Size!\n";
        throw(-1);
    }
    if (gamma.size() != compartmentDimensions[1])
    {
        Rcpp::Rcout << "Invalid gamma size!\n";
        throw(-1);
    }
    if (p_rs.size() != compartmentDimensions[1])
    {
        Rcpp::Rcout << "Invalid gamma size!\n";
        throw(-1);
    }



    Rcpp::Rcout << "Rcpp Provided Num Locations: " << compartmentDimensions[0] 
        << "\n";
    Rcpp::Rcout << "Rcpp Provided Num Times: " << compartmentDimensions[1] 
        << "\n";

    // Gather information for the creation of the 
    // covariate matrix
    covariateArgs xArgs;
    xArgs.inData_x = X.begin();
    xArgs.inData_z = Z.begin();
    xArgs.inRow_x = &covariateDimensions_x[0];
    xArgs.inCol_x = &covariateDimensions_x[1];
    xArgs.inRow_z = &covariateDimensions_z[0];
    xArgs.inCol_z = &covariateDimensions_z[1];

    // Gather information for the creation of the compartment matrices 
    
    compartmentArgs S_starArgs, E_starArgs, I_starArgs, R_starArgs;
    gammaArgs gammaFCArgs;
    S_starArgs.inData = S_star.begin();
    S_starArgs.inRow = &compartmentDimensions[0];
    S_starArgs.inCol = &compartmentDimensions[1];

    E_starArgs.inData = E_star.begin();
    E_starArgs.inRow = &compartmentDimensions[0];
    E_starArgs.inCol = &compartmentDimensions[1];

    I_starArgs.inData = I_star.begin();
    I_starArgs.inRow = &compartmentDimensions[0];
    I_starArgs.inCol = &compartmentDimensions[1];

    R_starArgs.inData = R_star.begin();
    R_starArgs.inRow = &compartmentDimensions[0];
    R_starArgs.inCol = &compartmentDimensions[1];

    gammaFCArgs.gamma = gamma.begin();
    gammaFCArgs.priorAlpha = priorAlpha_gamma.begin();
    gammaFCArgs.priorBeta = priorBeta_gamma.begin();

    // Gather information for the creation of the distance matrices

    double phi = 60*60*2.0;
    distanceArgs rawDistArgs; scaledDistanceArgs scaledDistArgs;
    rawDistArgs.inData = DistMat.begin(); 
    rawDistArgs.dim = &compartmentDimensions[0];
    scaledDistArgs.phi = &phi; 
    scaledDistArgs.inData = DistMat.begin();
    scaledDistArgs.dim = &compartmentDimensions[0];
    Rcpp::Rcout << "Loading covariate information into model context object\n";

    // Create the InitData object 
    InitData A0;
    A0.populate(S0.begin(),E0.begin(),I0.begin(),R0.begin(),
                S_star0.begin(),E_star0.begin(),I_star0.begin(),
                R_star0.begin(),&compartmentDimensions[0]);

    Rcpp::Rcout << "Creating Model Context\n";
    // Create the empty ModelContext object  
    ModelContext* context = new ModelContext();

    Rcpp::Rcout << "Populating Model Context\n";
    //Rcpp::Rcout << compartmentDimensions[0] << " " << compartmentDimensions[1] << "\n";
    //Rcpp::Rcout << (xArgs.inData_x)[1] << "\n";
    context -> populate(&A0, &xArgs, &S_starArgs, &E_starArgs, &I_starArgs, 
                        &R_starArgs, &rawDistArgs,&scaledDistArgs, &gammaFCArgs,
                        rho.begin(),beta.begin(),p_ei.begin(), p_ir.begin(),
                        p_rs.begin(),N.begin());

    // Set up output stream
    context -> fileProvider -> populate(context, chainOutputFile,
            (int*) chainOutputControl.begin(),(int*) chainStride.begin());
    context -> runSimulation_CPU(100000,*(verbose.begin()),*(debug.begin()));
    context -> fileProvider -> close();
    Rcpp::XPtr<ModelContext*> ptr(&context, true);
    // Clean up
    Rcpp::Rcout << "Finished.\n";
    delete chainOutputFile;
    return ptr;
}
