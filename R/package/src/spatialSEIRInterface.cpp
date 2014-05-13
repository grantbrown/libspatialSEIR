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

class spatialSEIRInterface
{

    private:
        //Attributes: 
        ModelContext* context;
        std::string* chainOutputFile;
        int* verbose;
        int* debug;

    public: 

        spatialSEIRInterface();
        int buildSpatialSEIRInterface(SEXP compMatDim,
                     SEXP xDim,
                     SEXP zDim,
                     SEXP xPrsDim,
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
                     SEXP X_pRS_,
                     SEXP DistMat_,
                     SEXP rho_,
                     SEXP gamma_,
                     SEXP priorAlpha_gamma_,
                     SEXP priorBeta_gamma_,
                     SEXP priorAlpha_pEI_,
                     SEXP priorBeta_pEI_,
                     SEXP priorAlpha_pIR_,
                     SEXP priorBeta_pIR_,
                     SEXP beta_,
                     SEXP betaPriorPrecision_,
                     SEXP betaPrs_,
                     SEXP betaPrsPriorPrecision_,
                     SEXP p_ei_,
                     SEXP p_ir_,
                     SEXP N_,
                     SEXP outFile,
                     SEXP logVarList,
                     SEXP iterationStride,
                     SEXP steadyStateConstraintPrecision_,
                     SEXP verboseFlag,
                     SEXP debugFlag,
                     SEXP sliceWidths);
        // Simulation Functions
        virtual int setRandomSeed(int seedVal);
        virtual int simulate(int iters);

        // Calculation Functions
        virtual int printDebugInfo();
        virtual int calculateS();
        virtual int calculateE();
        virtual int calculateI();
        virtual int calculateR();
        virtual int calculateP_SE(); 
        virtual int calculateP_RS();
        
        // Property Getter Functions
        // (All of this stuff is read only, should 
        // be changed only by calls to libspatialSEIR)
        
        virtual Rcpp::IntegerMatrix getS();
        virtual Rcpp::IntegerMatrix getE();
        virtual Rcpp::IntegerMatrix getI();
        virtual Rcpp::IntegerMatrix getR();

        virtual Rcpp::IntegerMatrix getS_star();
        virtual Rcpp::IntegerMatrix getE_star();
        virtual Rcpp::IntegerMatrix getI_star();
        virtual Rcpp::IntegerMatrix getR_star();

        virtual Rcpp::NumericMatrix getP_SE();
        virtual Rcpp::NumericVector getP_EI();
        virtual Rcpp::NumericVector getP_IR();
        virtual Rcpp::NumericVector getP_RS();
        virtual Rcpp::NumericVector getGamma();
        virtual Rcpp::NumericVector getBeta();
        virtual Rcpp::NumericVector getBetaP_RS();
        virtual Rcpp::NumericVector getRho();

        virtual int getDebug();
        virtual void setDebug(int debug_);

        virtual int getVerbose();
        virtual void setVerbose(int verbose_);

 
        //Destructor
        ~spatialSEIRInterface();
};

int spatialSEIRInterface::setRandomSeed(int seedVal)
{
    context -> setRandomSeed(seedVal);
    return(0);
}
int spatialSEIRInterface::simulate(int iters)
{
    context -> runSimulation_CPU(iters,*(verbose),*(debug));
    return(0);
}
int spatialSEIRInterface::printDebugInfo()
{
    Rcpp::Rcout << "S Dimensions: ";
    Rcpp::Rcout << *(context -> S -> nrow) << ", " << *(context -> S -> ncol) << "\n";
    return(0);
}

int spatialSEIRInterface::calculateS()
{
    if (*(context -> isPopulated))
    {
        (context -> calculateS_CPU());
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate S on non-populated ModelContext.\n";
    return(-1);
}
int spatialSEIRInterface::calculateE()
{
    if (*(context -> isPopulated))
    {
        (context -> calculateE_CPU());
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate E on non-populated ModelContext.\n";
    return(-1);
}
int spatialSEIRInterface::calculateI()
{
    if (*(context -> isPopulated))
    {
        (context -> calculateI_CPU());
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate I on non-populated ModelContext.\n";
    return(-1);
}
int spatialSEIRInterface::calculateR()
{
    if (*(context -> isPopulated))
    {
        (context -> calculateR_CPU());
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate R on non-populated ModelContext.\n";
    return(-1);
}
int spatialSEIRInterface::calculateP_SE() 
{
    if (*(context -> isPopulated))
    {
        (context -> calculateP_SE_CPU());
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate P_SE on non-populated ModelContext.\n";
    return(-1);
}
int spatialSEIRInterface::calculateP_RS()
{
    if (*(context -> isPopulated))
    {
        (context -> calculateP_RS_CPU());
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate P_RS on non-populated ModelContext.\n";
    return(-1);
}

Rcpp::IntegerMatrix spatialSEIRInterface::getS()
{
    Rcpp::IntegerMatrix output(*(context->S->nrow), *(context->S->ncol));
    int numVals = (*(context->S->nrow))*(*(context->S->ncol));
    int i;
    for (i = 0; i < numVals; i++)
    {
        output[i] = (context->S->data)[i];
    }
    return(output);
}
Rcpp::IntegerMatrix spatialSEIRInterface::getE()
{
    Rcpp::IntegerMatrix output(*(context->S->nrow), *(context->S->ncol));
    int numVals = (*(context->S->nrow))*(*(context->S->ncol));
    int i;
    for (i = 0; i < numVals; i++)
    {
        output[i] = (context->E->data)[i];
    }
    return(output);

}
Rcpp::IntegerMatrix spatialSEIRInterface::getI()
{
    Rcpp::IntegerMatrix output(*(context->S->nrow), *(context->S->ncol));
    int numVals = (*(context->S->nrow))*(*(context->S->ncol));
    int i;
    for (i = 0; i < numVals; i++)
    {
        output[i] = (context->I->data)[i];
    }
    return(output);

}
Rcpp::IntegerMatrix spatialSEIRInterface::getR()
{
    Rcpp::IntegerMatrix output(*(context->S->nrow), *(context->S->ncol));
    int numVals = (*(context->S->nrow))*(*(context->S->ncol));
    int i;
    for (i = 0; i < numVals; i++)
    {
        output[i] = (context->R->data)[i];
    }
    return(output);

}

Rcpp::IntegerMatrix spatialSEIRInterface::getS_star()
{
    Rcpp::IntegerMatrix output(*(context->S->nrow), *(context->S->ncol));
    int numVals = (*(context->S->nrow))*(*(context->S->ncol));
    int i;
    for (i = 0; i < numVals; i++)
    {
        output[i] = (context->S_star->data)[i];
    }
    return(output);

}
Rcpp::IntegerMatrix spatialSEIRInterface::getE_star()
{
    Rcpp::IntegerMatrix output(*(context->S->nrow), *(context->S->ncol));
    int numVals = (*(context->S->nrow))*(*(context->S->ncol));
    int i;
    for (i = 0; i < numVals; i++)
    {
        output[i] = (context->E_star->data)[i];
    }
    return(output);

}
Rcpp::IntegerMatrix spatialSEIRInterface::getI_star()
{
    Rcpp::IntegerMatrix output(*(context->S->nrow), *(context->S->ncol));
    int numVals = (*(context->S->nrow))*(*(context->S->ncol));
    int i;
    for (i = 0; i < numVals; i++)
    {
        output[i] = (context->I_star->data)[i];
    }
    return(output);

}
Rcpp::IntegerMatrix spatialSEIRInterface::getR_star()
{
    Rcpp::IntegerMatrix output(*(context->S->nrow), *(context->S->ncol));
    int numVals = (*(context->S->nrow))*(*(context->S->ncol));
    int i;
    for (i = 0; i < numVals; i++)
    {
        output[i] = (context->R_star->data)[i];
    }
    return(output);

}

Rcpp::NumericMatrix spatialSEIRInterface::getP_SE()
{
    Rcpp::NumericMatrix output(*(context->S->nrow), *(context->S->ncol));
    int numVals = (*(context->S->nrow))*(*(context->S->ncol));
    int i;
    for (i = 0; i < numVals; i++)
    {
        output[i] = (context->p_se)[i];
    }
    return(output);

}
Rcpp::NumericVector spatialSEIRInterface::getP_RS()
{
    Rcpp::NumericVector output(*(context->S->ncol));
    int i;
    int numVals = (*(context->S->ncol));
    for (i = 0; i < numVals; i++) 
    {
        output[i] = (context->p_rs)[i]; 
    }
    return(output);
}
Rcpp::NumericVector spatialSEIRInterface::getGamma()
{
    Rcpp::NumericVector output(*(context->S->ncol));
    int i;
    int numVals = (*(context->S->ncol));
    for (i = 0; i < numVals; i++) 
    {
        output[i] = (context->gamma)[i]; 
    }
    return(output);
}
Rcpp::NumericVector spatialSEIRInterface::getBeta()
{
    Rcpp::NumericVector output(((*(context->X->ncol_x)) + (*(context->X->ncol_z))));
    int i;
    int numVals = (((*(context->X->ncol_x)) + (*(context->X->ncol_z))));
    for (i = 0; i < numVals; i++) 
    {
        output[i] = (context->beta)[i]; 
    }
    return(output);

}
Rcpp::NumericVector spatialSEIRInterface::getBetaP_RS()
{
    Rcpp::NumericVector output(*(context->X_pRS->ncol_x));
    int i;
    int numVals = (*(context->X_pRS->ncol_x));
    for (i = 0; i < numVals; i++) 
    {
        output[i] = (context->betaPrs)[i]; 
    }
    return(output);
}
Rcpp::NumericVector spatialSEIRInterface::getP_EI()
{
    Rcpp::NumericVector output(1);
    output[0] = *(context->p_ei); 
    return(output);
}


Rcpp::NumericVector spatialSEIRInterface::getP_IR()
{
    Rcpp::NumericVector output(1);
    output[0] = *(context->p_ir); 
    return(output);
}

Rcpp::NumericVector spatialSEIRInterface::getRho()
{
    Rcpp::NumericVector output(1);
    output[0] = *(context->rho); 
    return(output);
}

int spatialSEIRInterface::getDebug()
{
    int out;
    out = *debug;
    return(out);    
}
void spatialSEIRInterface::setDebug(int debug_)
{
   *debug = debug_;  
}
int spatialSEIRInterface::getVerbose()
{
    int output;
    output= *verbose;
    return(output);
}
void spatialSEIRInterface::setVerbose(int verbose_)
{
   *verbose = verbose_; 
}

spatialSEIRInterface::spatialSEIRInterface()
{
    Rcpp::Rcout << "Creating Model Context\n";
    // Create the empty ModelContext object  
    context = new ModelContext();
    verbose = new int();
    debug = new int(); 
    *verbose = 0;
    *debug = 0;
}

int spatialSEIRInterface::buildSpatialSEIRInterface(SEXP compMatDim,
                     SEXP xDim,
                     SEXP zDim,
                     SEXP xPrsDim,
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
                     SEXP X_pRS_,
                     SEXP DistMat_,
                     SEXP rho_,
                     SEXP gamma_,
                     SEXP priorAlpha_gamma_,
                     SEXP priorBeta_gamma_,
                     SEXP priorAlpha_pEI_,
                     SEXP priorBeta_pEI_,
                     SEXP priorAlpha_pIR_,
                     SEXP priorBeta_pIR_,
                     SEXP beta_,
                     SEXP betaPriorPrecision_,
                     SEXP betaPrs_,
                     SEXP betaPrsPriorPrecision_,
                     SEXP p_ei_,
                     SEXP p_ir_,
                     SEXP N_,
                     SEXP outFile,
                     SEXP logVarList,
                     SEXP iterationStride,
                     SEXP steadyStateConstraintPrecision_,
                     SEXP verboseFlag,
                     SEXP debugFlag,
                     SEXP sliceWidths)
{
    int err = 0;
    Rcpp::Rcout << "Wrapping input data in Rcpp vectors.\n";
    //Deal with the data conversion from R to c++
    Rcpp::IntegerVector compartmentDimensions(compMatDim);
    Rcpp::IntegerVector covariateDimensions_x(xDim);
    Rcpp::IntegerVector covariateDimensions_z(zDim);
    Rcpp::IntegerVector covariateDimension_pRS_x(xPrsDim);
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
    Rcpp::NumericVector X_pRS(X_pRS_);
    Rcpp::NumericVector DistMat(DistMat_);

    Rcpp::NumericVector rho(rho_);

    Rcpp::NumericVector gamma(gamma_);
    Rcpp::NumericVector priorAlpha_gamma(priorAlpha_gamma_);
    Rcpp::NumericVector priorBeta_gamma(priorBeta_gamma_);

    Rcpp::NumericVector priorAlpha_pEI(priorAlpha_pEI_);
    Rcpp::NumericVector priorBeta_pEI(priorBeta_pEI_);
    Rcpp::NumericVector priorAlpha_pIR(priorAlpha_pIR_);
    Rcpp::NumericVector priorBeta_pIR(priorBeta_pIR_);


    Rcpp::NumericVector beta(beta_);
    Rcpp::NumericVector betaPriorPrecision(betaPriorPrecision_);
    Rcpp::NumericVector betaPrs(betaPrs_);
    Rcpp::NumericVector betaPrsPriorPrecision(betaPrsPriorPrecision_);
    Rcpp::NumericVector p_ei(p_ei_);
    Rcpp::NumericVector p_ir(p_ir_);
    Rcpp::IntegerVector N(N_);

    Rcpp::NumericVector steadyStateConstraintPrecision(steadyStateConstraintPrecision_);
    Rcpp::NumericVector sliceParams(sliceWidths);

    chainOutputFile = new std::string(); 
    *chainOutputFile = Rcpp::as<std::string>(outFile);
    Rcpp::IntegerVector chainOutputControl(logVarList); 

    Rcpp::IntegerVector vFlag(verboseFlag);
    Rcpp::IntegerVector dFlag(debugFlag);
    *verbose = vFlag[0];
    *debug = dFlag[0];
    // logVarList: (Beta, rho,gamma,p_se,p_ei,p_ir,betaPrs,S*, E*, I*, R*)
    // Nonzero if respective variables are to be output
    Rcpp::IntegerVector chainStride(iterationStride);
    

    // Sanity check the input data. 
    if (compartmentDimensions.size() != 2)
    {
        Rcpp::Rcout << "Compartments must be two dimensional.\n";
        throw(-1);
    }
    if (covariateDimensions_x.size() != 2 || covariateDimensions_z.size() != 2)
    {
        Rcpp::Rcout << "Covariates must be two dimensional.\n";
        throw(-1);
    }
    if (vFlag.size() != 1 || dFlag.size() != 1)
    {
        Rcpp::Rcout << "Verbose and debug flags must be length 1\n";
    }
    if (chainStride.size() != 1)
    {
        Rcpp::Rcout << "Chain stride must be length 1\n";
    }

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
    if (N.size() != compartmentSize)
    {
        Rcpp::Rcout << "Invalid N Compartment Size!\n";
        throw(-1);
    }
    if ((X_pRS.size() % compartmentDimensions[1]) != 0)
    {
        Rcpp::Rcout << "Invalid X_pRS size.\n";
        Rcpp::Rcout << "Size: " << X_pRS.size() << ", Number of Time Points: " << compartmentDimensions[1] << "\n";
    }

    if (gamma.size() != compartmentDimensions[1])
    {
        Rcpp::Rcout << "Invalid gamma size!\n";
        Rcpp::Rcout << "Size: " << gamma.size() << ", Number of Time Points: " << compartmentDimensions[1] << "\n";
        throw(-1);
    }
    if (sliceParams.size() != 7)
    {
        Rcpp::Rcout << "Slice sampling parameters must be of length 6: S*,E*,R*,beta,betaPrs,rho,gamma\n";
        throw(-1);
    }
    if (chainOutputControl.size() != 30)
    {
        Rcpp::Rcout << "There are 30 chain output options, please don't leave any blank.\n";
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

    covariateArgs xPrsArgs; 
    xPrsArgs.inData_x = X_pRS.begin();
    xPrsArgs.inData_z = NULL;
    xPrsArgs.inRow_x = &covariateDimension_pRS_x[0];
    xPrsArgs.inCol_x = &covariateDimension_pRS_x[1];
    // Clean this up, pass values instead. 
    int zeroVal = 0;
    xPrsArgs.inRow_z = &zeroVal;
    xPrsArgs.inCol_z = &zeroVal;


    // Gather information for the creation of the compartment matrices 
    
    compartmentArgs S_starArgs, E_starArgs, I_starArgs, R_starArgs;
    gammaArgs gammaFCArgs;
    sliceParameters sliceParamStruct;

    sliceParamStruct.S_starWidth = &sliceParams[0];
    sliceParamStruct.E_starWidth = &sliceParams[1];
    sliceParamStruct.R_starWidth = &sliceParams[2];
    sliceParamStruct.betaWidth = &sliceParams[3];
    sliceParamStruct.betaPrsWidth = &sliceParams[4];
    sliceParamStruct.rhoWidth = &sliceParams[5];
    sliceParamStruct.gammaWidth = &sliceParams[6];

    S_starArgs.inData = S_star.begin();
    S_starArgs.inRow = &compartmentDimensions[0];
    S_starArgs.inCol = &compartmentDimensions[1];
    S_starArgs.steadyStateConstraintPrecision = steadyStateConstraintPrecision[0];

    E_starArgs.inData = E_star.begin();
    E_starArgs.inRow = &compartmentDimensions[0];
    E_starArgs.inCol = &compartmentDimensions[1];
    E_starArgs.steadyStateConstraintPrecision = steadyStateConstraintPrecision[0];


    I_starArgs.inData = I_star.begin();
    I_starArgs.inRow = &compartmentDimensions[0];
    I_starArgs.inCol = &compartmentDimensions[1];

    R_starArgs.inData = R_star.begin();
    R_starArgs.inRow = &compartmentDimensions[0];
    R_starArgs.inCol = &compartmentDimensions[1];
    R_starArgs.steadyStateConstraintPrecision = steadyStateConstraintPrecision[0];


    gammaFCArgs.gamma = gamma.begin();
    gammaFCArgs.priorAlpha = priorAlpha_gamma.begin();
    gammaFCArgs.priorBeta = priorBeta_gamma.begin();

    priorControl priorValues;
    priorValues.betaPriorPrecision = betaPriorPrecision[0];
    priorValues.P_EI_priorAlpha = priorAlpha_pEI[0];
    priorValues.P_EI_priorBeta = priorBeta_pEI[0];
    priorValues.P_IR_priorAlpha = priorAlpha_pIR[0];
    priorValues.P_IR_priorBeta = priorBeta_pIR[0];
    priorValues.betaPrsPriorPrecision = betaPrsPriorPrecision[0];

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

    Rcpp::Rcout << "Populating Model Context\n";
    //Rcpp::Rcout << compartmentDimensions[0] << " " << compartmentDimensions[1] << "\n";
    //Rcpp::Rcout << (xArgs.inData_x)[1] << "\n";
    context -> populate(&A0, &xArgs, &xPrsArgs, &S_starArgs, &E_starArgs, &I_starArgs, 
                        &R_starArgs, &rawDistArgs,&scaledDistArgs, &gammaFCArgs,
                        rho.begin(),beta.begin(),p_ei.begin(), p_ir.begin(),
                        betaPrs.begin(),N.begin(),&sliceParamStruct, &priorValues);

    // Set up output stream
    context -> fileProvider -> populate(context, chainOutputFile,
            (int*) chainOutputControl.begin(),(int*) chainStride.begin());

    Rcpp::Rcout << "Context setup complete.\n";
    return(err);
}


spatialSEIRInterface::~spatialSEIRInterface()
{   
    // Context handles the complicated cleanup
    delete chainOutputFile;
    delete context;
    delete verbose;
    delete debug;
}


RCPP_MODULE(mod_spatialSEIRInterface)
{
    using namespace Rcpp;
    class_<spatialSEIRInterface>( "spatialSEIRInterface" )

    .constructor()

    .method("buildSpatialSEIRInterface", &spatialSEIRInterface::buildSpatialSEIRInterface)
    .method("printDebugInfo", &spatialSEIRInterface::printDebugInfo)
    .method("setRandomSeed", &spatialSEIRInterface::setRandomSeed)
    .method("simulate", &spatialSEIRInterface::simulate)
    .method("calculateS", &spatialSEIRInterface::calculateS)
    .method("calculateE", &spatialSEIRInterface::calculateE)
    .method("calculateI", &spatialSEIRInterface::calculateI)
    .method("calculateR", &spatialSEIRInterface::calculateR)
    .method("calculateP_RS", &spatialSEIRInterface::calculateP_RS)
    .method("calculateP_SE", &spatialSEIRInterface::calculateP_SE)
    .property("S", &spatialSEIRInterface::getS, "Susceptible Compartment Matrix")
    .property("E", &spatialSEIRInterface::getE, "Exposed Compartment Matrix")
    .property("I", &spatialSEIRInterface::getI, "Infectious Compartment Matrix")
    .property("R", &spatialSEIRInterface::getR, "Removed Compartment Matrix")
    .property("S_star", &spatialSEIRInterface::getS_star, "Removed to Susceptible Transition Matrix")
    .property("E_star", &spatialSEIRInterface::getE_star, "Susceptible to Exposed Transition Matrix")
    .property("I_star", &spatialSEIRInterface::getI_star, "Exposed to Infectious Transition Matrix")
    .property("R_star", &spatialSEIRInterface::getR_star, "Infectious to Removed Transition Matrix")
    .property("p_se", &spatialSEIRInterface::getP_SE, "Exposure Probability Matrix")
    .property("p_ei", &spatialSEIRInterface::getP_EI, "E to I Transition Probability")
    .property("p_ir", &spatialSEIRInterface::getP_IR, "I to R Transition Probability")
    .property("p_rs", &spatialSEIRInterface::getP_RS, "R-S Transition Probability Vector")
    .property("gamma", &spatialSEIRInterface::getGamma, "External Infection Component Vector")
    .property("beta", &spatialSEIRInterface::getBeta, "Exposure Process Regression Parameters")
    .property("betaP_RS", &spatialSEIRInterface::getBetaP_RS, "R-S Transition Process Regression Parameters")
    .property("rho", &spatialSEIRInterface::getRho, "Spatial Dependence Term")
    .property("debug", &spatialSEIRInterface::getDebug, &spatialSEIRInterface::setDebug, "Show debug level output?")
    .property("verbose", &spatialSEIRInterface::getVerbose, &spatialSEIRInterface::setVerbose, "Sho verbose level output?")
    ;
}

