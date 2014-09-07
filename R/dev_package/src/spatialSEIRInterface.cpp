#include <Rcpp.h>
#include <cmath>
#include <ModelContext.hpp>
#include <LSS_FullConditionalList.hpp>
#include <LSS_Samplers.hpp>
#include <CovariateMatrix.hpp>
#include <CompartmentalModelMatrix.hpp>
#include <DistanceMatrix.hpp>
#include <IOProvider.hpp>
#include <OCLProvider.hpp>

using namespace Rcpp;
using namespace SpatialSEIR;

struct LSSCout {};
extern LSSCout lssCout;

template <typename T>
    LSSCout& operator<< (LSSCout &s, const T &x){
        Rcpp::Rcout << x;
        return s;
    }

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
                     SEXP Y,
                     SEXP Sstar, 
                     SEXP Estar, 
                     SEXP Istar, 
                     SEXP Rstar, 
                     SEXP offset_,
                     SEXP X_,
                     SEXP Z_,
                     SEXP X_pRS_,
                     SEXP DistMat_,
                     SEXP rho_,
                     SEXP priorAlpha_pEI_,
                     SEXP priorBeta_pEI_,
                     SEXP priorAlpha_pIR_,
                     SEXP priorBeta_pIR_,
                     SEXP beta_,
                     SEXP betaPriorPrecision_,
                     SEXP betaPrs_,
                     SEXP betaPrsPriorPrecision_,
                     SEXP gamma_ei_,
                     SEXP gamma_ir_,
                     SEXP N_,
                     SEXP outFile,
                     SEXP iterationStride,
                     SEXP steadyStateConstraintPrecision_,
                     SEXP verboseFlag,
                     SEXP debugFlag,
                     SEXP sliceWidths,
                     SEXP reinfectionMode,
                     SEXP scaleDistanceMode_);
        // Simulation Functions
        virtual int setRandomSeed(int seedVal);
        virtual int simulate(int iters);
        virtual void setPredictionTraces();
        virtual int setTrace(int locationIndex);
        virtual int setTrace2(int locationIndex, int timeIndex);
        virtual void setDevice(int platformId, int deviceId);
        virtual void setCompartmentSamplingMode(int mode);
        virtual int getCompartmentSamplingMode();
        virtual void setParameterSamplingMode(int mode);
        virtual int getParameterSamplingMode();



        // Calculation Functions
        virtual int printDebugInfo();
        virtual int calculateS();
        virtual int calculateE();
        virtual int calculateI();
        virtual int calculateR();
        virtual int calculateP_SE(); 
        virtual int calculateP_SE2(int i, int j); 
        virtual int calculateP_SE_OCL();
        virtual double estimateR0();
        virtual double estimateR02(int t);
        virtual double estimateR03(int i, int t);
        virtual int calculateP_RS();

        
        // Property Getter Functions
        // (All of this stuff is read only, should 
        // be changed only by calls to libspatialSEIR)
        virtual void updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange);
        virtual void printSamplingParameters();
        virtual void printAcceptanceRates();       
        virtual void printOCLSummary();
        virtual Rcpp::IntegerMatrix getS();
        virtual Rcpp::IntegerMatrix getE();
        virtual Rcpp::IntegerMatrix getI();
        virtual Rcpp::IntegerMatrix getR();

        virtual Rcpp::IntegerMatrix getS_star();
        virtual Rcpp::IntegerMatrix getE_star();
        virtual Rcpp::IntegerMatrix getI_star();
        virtual Rcpp::IntegerMatrix getR_star();

        virtual Rcpp::IntegerVector getS0();
        virtual Rcpp::IntegerVector getE0();
        virtual Rcpp::IntegerVector getI0();
        virtual Rcpp::IntegerVector getR0();

        virtual Rcpp::NumericMatrix getP_SE();
        virtual Rcpp::NumericVector getP_EI();
        virtual Rcpp::NumericVector getP_IR();
        virtual Rcpp::NumericVector getP_RS();
        virtual Rcpp::NumericVector getGenerationMatrix(int t);
        virtual Rcpp::NumericVector getIntegratedGenerationMatrix(int t);
        virtual Rcpp::NumericVector getBeta();
        virtual Rcpp::NumericVector getBetaP_RS();
        virtual Rcpp::NumericVector getRho();

        virtual int getDebug();
        virtual void setDebug(int debug_);

        virtual int getVerbose();
        virtual void setVerbose(int verbose_);

        virtual int getUseDecorrelation();
        virtual void setUseDecorrelation(int val);

        virtual void standardizeDistanceMatrix();
 
        //Destructor
        ~spatialSEIRInterface();
};


int spatialSEIRInterface::getUseDecorrelation()
{
    if (*(context -> isPopulated))
    {
        return((context -> config -> useDecorrelation));
    }
    Rcpp::Rcout << "Context Not populated\n";
    return(-1);
}

void spatialSEIRInterface::setUseDecorrelation(int val)
{
    if (*(context -> isPopulated))
    {
        ((context -> config ->useDecorrelation) = val);
        (context -> configureIterationTasks());
        return;
    }
    Rcpp::Rcout << "Context Not populated\n";
    return; 
}

int spatialSEIRInterface::setRandomSeed(int seedVal)
{
    if (*(context -> isPopulated))
    {
        context -> setRandomSeed(seedVal);
        return(0);
    }
    Rcpp::Rcout << "Context Not populated\n";
    return(-1);
}
int spatialSEIRInterface::simulate(int iters)
{
    if (*(context -> isPopulated))
    {
        context -> runSimulation(iters,*(verbose),*(debug));
        return(0);
    }
    Rcpp::Rcout << "Context Not populated\n";

}

void spatialSEIRInterface::setDevice(int platformId, int deviceId)
{
    if (*(context -> isPopulated))
    {
        (context -> oclProvider -> setDevice(platformId, deviceId));
        return;
    }
    Rcpp::Rcout << "ModelContext is not populated.\n";
    return;
}

void spatialSEIRInterface::setCompartmentSamplingMode(int mode)
{
    int oldMode = getCompartmentSamplingMode();
    if (*(context -> isPopulated))
    {
        if (mode == COMPARTMENT_METROPOLIS_SAMPLER)
        {
            Rcpp::Rcout << "Setting compartment sampling mode to Metropolis\n";
        }
        else if (mode == COMPARTMENT_IDX_METROPOLIS_SAMPLER)
        {
            Rcpp::Rcout << "Setting compartment sampling mode to indexed Metropolis\n";
        }
        else if (mode == COMPARTMENT_IDX_SLICE_SAMPLER)
        {
            Rcpp::Rcout << "Setting compartment sampling mode to indexed slice\n";
        }
        else if (mode == COMPARTMENT_METROPOLIS_SAMPLER_OCL)
        {
             Rcpp::Rcout << "Setting compartment sampling mode to Metropolis using OpenCL\n";  
        }
        else if (mode == COMPARTMENT_BINOM_PROPOSAL_METROPOLIS_SAMPLER)
        {
             Rcpp::Rcout << "Setting compartment sampling mode to Metropolis with chain binomial proposal.\n";  
        }
        else if (mode == COMPARTMENT_BINOM_PROPOSAL_SLICE_SAMPLER)
        {
             Rcpp::Rcout << "Setting compartment sampling mode to slice with chain binomial proposal.\n";  
        }
        else
        {
            Rcpp::Rcout << "Error: mode must be one of the following: \n"
                        << COMPARTMENT_METROPOLIS_SAMPLER << ": Metropolis \n   "
                        << COMPARTMENT_METROPOLIS_SAMPLER_OCL << ": Metropolis with OpenCL\n   "
                        << COMPARTMENT_IDX_METROPOLIS_SAMPLER << ": Indexed Metropolis\n   "
                        << COMPARTMENT_BINOM_PROPOSAL_METROPOLIS_SAMPLER << ": Metropolis w/ chain binomial proposal.\n" 
                        << COMPARTMENT_BINOM_PROPOSAL_SLICE_SAMPLER << ": Slice w/ chain binomial proposal.\n" 
                        << COMPARTMENT_IDX_SLICE_SAMPLER << ": Indexed Slice\n   ";
            return;
        }

        try
        {
            context -> setCompartmentSamplingMode(mode);
        }
        catch (int err)
        {
            Rcpp::Rcout << "Unable to update compartment sampling mode\n";
            setCompartmentSamplingMode(oldMode);
        }
        return;
    }
    Rcpp::Rcout << "Context Not populated\n";
}

int spatialSEIRInterface::getCompartmentSamplingMode()
{
    if (*(context -> isPopulated))
    {
        return((context -> getCompartmentSamplingMode()));
    }
    return(-1);
}

void spatialSEIRInterface::setParameterSamplingMode(int mode)
{

    int oldMode = getCompartmentSamplingMode();
    if (*(context -> isPopulated))
    {
        if (mode == PARAMETER_SINGLE_METROPOLIS_SAMPLER)
        {
            Rcpp::Rcout << "Setting parameter sampling mode to single Metropolis\n";
        }
        else if (mode == PARAMETER_JOINT_METROPOLIS_SAMPLER)
        {
            Rcpp::Rcout << "Setting parameter sampling mode to joint Metropolis\n";
        }
        else if (mode == PARAMETER_JOINT_SLICE_SAMPLER)
        {
            Rcpp::Rcout << "Setting parameter sampling mode to joint slice\n";
        }
        else if (mode == PARAMETER_DECORR_SAMPLER)
        {
            Rcpp::Rcout << "Setting parameter sampling mode to decorrelation Metropolis\n";
        }
        else if (mode == PARAMETER_JOINT_METROPOLIS_SAMPLER_OCL)
        {
            Rcpp::Rcout << "Setting parameter sampling mode to joint Metropolis with OpenCL\n";
        }
        else
        {
            Rcpp::Rcout << "Error: mode must be one of the following: \n"
                        << PARAMETER_SINGLE_METROPOLIS_SAMPLER << ": Single Value Metropolis \n   "
                        << PARAMETER_JOINT_METROPOLIS_SAMPLER << ": Joint Proposal Metropolis\n   "
                        << PARAMETER_JOINT_SLICE_SAMPLER << ": Joint Slice\n   "
                        << PARAMETER_DECORR_SAMPLER << ": Decorrelation Proposal Metropolis\n   "
                        << PARAMETER_JOINT_METROPOLIS_SAMPLER_OCL << ": Joint Proposal Metropolis with OpenCL\n   "                        
                        ;
            return;
        }

        try
        {
            context -> setParameterSamplingMode(mode);
        }
        catch (int err)
        {
            Rcpp::Rcout << "Unable to update compartment sampling mode\n";
            setParameterSamplingMode(oldMode);
        }
        return;
    }
    Rcpp::Rcout << "Context Not populated\n";
}

int spatialSEIRInterface::getParameterSamplingMode()
{
    if (*(context -> isPopulated))
    {
        return((context -> getParameterSamplingMode()));
    }
    return(-1);
}



void spatialSEIRInterface::standardizeDistanceMatrix()
{
    if (*(context -> isPopulated))
    {
       if (*(context -> numIterations) != 0)
       {
           Rcpp::Rcout << "Can't change distance matrix once sampling has begun.\n";
           return;
       }
       (context -> scaledDistMat -> makeRowStochastic()); 
    }
    else
    {
        Rcpp::Rcout << "No distance matrix to standardize.\n";
    }
}


int spatialSEIRInterface::setTrace(int locationIndex)
{
    if (*(context -> isPopulated))
    {
        try
        {
            (context -> fileProvider -> setTrace(locationIndex));
        }
        catch (int err)
        {
            Rcpp::Rcout << "Unable to set trace for location: " << locationIndex << "\n";
        }
        return(0);
    }
    Rcpp::Rcout << "Attept to set trace on non-populated ModelContext.\n";
    return(-1);
}
int spatialSEIRInterface::setTrace2(int locationIndex, int timeIndex)
{
    if (*(context -> isPopulated))
    {
        try
        {
            (context -> fileProvider -> setTrace(locationIndex, timeIndex));
        }
        catch (int err)
        {
            Rcpp::Rcout << "Unable to set trace for (loc,time): (" << locationIndex 
                        << ", " << timeIndex << ")"<< "\n";
        }

        return(0);
    }
    Rcpp::Rcout << "Attept to set trace on non-populated ModelContext.\n";
    return(-1);
}

void spatialSEIRInterface::setPredictionTraces()
{
    if (*(context -> isPopulated))
    {
        int i;
        int lastTimeIdx = (*(context -> S_star -> nrow))-1;
        for (i = 0; i < *(context -> S_star -> ncol); i++)
        {
            setTrace2(i,lastTimeIdx);
        }
        return;
    }
    Rcpp::Rcout << "Attept to set trace on non-populated ModelContext.\n";
    return;

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
int spatialSEIRInterface::calculateP_SE2(int i, int j) 
{
    if (*(context -> isPopulated))
    {
        (context -> calculateP_SE_CPU(i,j));
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate P_SE on non-populated ModelContext.\n";
    return(-1);
}
int spatialSEIRInterface::calculateP_SE_OCL() 
{
    if (*(context -> isPopulated))
    {
        (context -> calculateP_SE_OCL());
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate P_SE on non-populated ModelContext.\n";
    return(-1);
}

double spatialSEIRInterface::estimateR0()
{
    if (*(context -> isPopulated))
    {
        return((context -> estimateR0()));
    }
    Rcpp::Rcout << "Attempt to estimate R0 on a non-populated ModelContext.\n";
    return(-1.0);
}

double spatialSEIRInterface::estimateR02(int t)
{
    if (*(context -> isPopulated))
    {
        return((context -> estimateR0(t)));
    }
    Rcpp::Rcout << "Attempt to estimate R0 on a non-populated ModelContext.\n";
    return(-1.0);
}

double spatialSEIRInterface::estimateR03(int i, int t)
{
    if (*(context -> isPopulated))
    {
        return((context -> estimateR0(i, t)));
    }
    Rcpp::Rcout << "Attempt to estimate R0 on a non-populated ModelContext.\n";
    return(-1.0);
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

void spatialSEIRInterface::updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange)
{
    if ((*(context -> numIterations)) == 0)
    {
        Rcpp::Rcout << "No samples drawn.\n";
        return;
    }
    (context -> updateSamplingParameters(desiredRatio, targetWidth, proportionChange));
}

void spatialSEIRInterface::printOCLSummary()
{
    if (*(context -> isPopulated))
    {
        (context -> oclProvider -> printSummary());
        return;
    }
    Rcpp::Rcout << "ModelContext has not been populated.\n";
}

void spatialSEIRInterface::printSamplingParameters()
{
    int i;
    Rcpp::Rcout << "S0:       " << (*(context -> S0_fc -> sliceWidth)) << "\n"; 
    Rcpp::Rcout << "I0:       " << (*(context -> I0_fc -> sliceWidth)) << "\n"; 
    if ((context -> config -> reinfectionMode) <= 2)
    {
        Rcpp::Rcout << "S_star:   " << (*(context -> S_star_fc -> sliceWidth)) << "\n"; 
    }
    Rcpp::Rcout << "E_star:   " << (*(context -> E_star_fc -> sliceWidth)) << "\n"; 
    Rcpp::Rcout << "R_star:   " << (*(context -> R_star_fc -> sliceWidth))<< "\n";  

    Rcpp::Rcout << "beta:     ";
    for (i = 0; i < (*(context -> beta_fc -> varLen)); i++)
    {
        Rcpp::Rcout << ((context -> beta_fc -> sliceWidth)[i]) << ", ";

    }
    Rcpp::Rcout << "\n"; 

    if ((context -> config -> reinfectionMode) == 1)
    {
        Rcpp::Rcout << "betaP_RS:     ";
        for (i = 0; i < (*(context -> betaPrs_fc -> varLen)); i++)
        {
            Rcpp::Rcout << ((context -> betaPrs_fc -> sliceWidth)[i]) << ", ";

        }
        Rcpp::Rcout << "\n"; 
    }

    if (*(context -> S_star -> ncol) > 1)
    {
        Rcpp::Rcout << "rho:      " << (*(context -> rho_fc -> sliceWidth)) << "\n"; 
    }
    Rcpp::Rcout << "gamma_ei: " << (*(context -> gamma_ei_fc -> sliceWidth)) << "\n"; 
    Rcpp::Rcout << "gamma_ir: " <<  (*(context -> gamma_ir_fc -> sliceWidth)) << "\n"; 

}

void spatialSEIRInterface::printAcceptanceRates()
{
    int i;
    if ((*(context -> numIterations)) == 0)
    {
        Rcpp::Rcout << "No samples drawn.\n";
        return;
    }
    Rcpp::Rcout << "Total iterations so far: " << (*(context -> numIterations)) << "\n";
    Rcpp::Rcout << "Acceptance rates: \n";
    Rcpp::Rcout << "S0:       " << (*(context -> S0_fc -> accepted)*1.0)/(*(context -> S0_fc -> samples)) 
                          << "\n"; 
    Rcpp::Rcout << "I0:       " << (*(context -> I0_fc -> accepted)*1.0)/ 
                                  (*(context -> I0_fc -> samples)) 
                          <<  "\n"; 
    if ((context -> config -> reinfectionMode) <= 2)
    {
        Rcpp::Rcout << "S_star:   " << (*(context -> S_star_fc -> accepted)*1.0)/
                                          (*(context -> S_star_fc -> samples)) 
                                  << "\n"; 
    }

    Rcpp::Rcout << "E_star:   " << (*(context -> E_star_fc -> accepted)*1.0)/
                                      (*(context -> E_star_fc -> samples)) 
                              << "\n"; 
    Rcpp::Rcout << "R_star:   " << (*(context -> R_star_fc -> accepted)*1.0)/
                                      (*(context -> R_star_fc -> samples)) 
                              << "\n";  
    Rcpp::Rcout << "beta:     ";
    for (i = 0; i < *(context -> beta_fc -> varLen); i++)
    {
        Rcpp::Rcout << ((context -> beta_fc -> accepted)[i]*1.0)/
                                      (*(context -> beta_fc -> samples))
                                      << ", "; 
    }
    Rcpp::Rcout << "\n"; 

    if (*(context -> S_star -> ncol) > 1)
    {
        Rcpp::Rcout << "rho:      " << (*(context -> rho_fc -> accepted)*1.0)/
                                          (*(context -> rho_fc -> samples)) 
                                  << "\n"; 
    }
    if ((context -> config -> reinfectionMode) == 1)
    {
        Rcpp::Rcout << "betaP_RS: ";
        for (i = 0; i < *(context -> betaPrs_fc -> varLen); i++)
        {
            Rcpp::Rcout << ((context -> betaPrs_fc -> accepted)[i]*1.0)/
                                      (*(context -> betaPrs_fc -> samples)) 
                                      << ", "; 
        }
        Rcpp::Rcout << "\n"; 
    }
    Rcpp::Rcout << "gamma_ei:     " << (*(context -> gamma_ei_fc -> accepted)*1.0)/
                                      (*(context -> gamma_ei_fc -> samples)) 
                              << "\n"; 

    Rcpp::Rcout << "gamma_ir:     " << (*(context -> gamma_ir_fc -> accepted)*1.0)/
                                      (*(context -> gamma_ir_fc -> samples))
                              << "\n"; 

}

Rcpp::IntegerVector spatialSEIRInterface::getS0()
{
    Rcpp::IntegerVector output(*(context->S->ncol));
    int numVals = *(context->S->ncol);
    int i;
    for (i = 0; i < numVals; i++)
    {
        output[i] = (context->A0 -> S0)[i];
    }
    return(output);
}
Rcpp::IntegerVector spatialSEIRInterface::getE0()
{
    Rcpp::IntegerVector output(*(context->S->ncol));
    int numVals = *(context->S->ncol);
    int i;
    for (i = 0; i < numVals; i++)
    {
        output[i] = (context->A0 -> E0)[i];
    }
    return(output);
}
Rcpp::IntegerVector spatialSEIRInterface::getI0()
{
    Rcpp::IntegerVector output(*(context->S->ncol));
    int numVals = *(context->S->ncol);
    int i;
    for (i = 0; i < numVals; i++)
    {
        output[i] = (context->A0 -> I0)[i];
    }
    return(output);
}
Rcpp::IntegerVector spatialSEIRInterface::getR0()
{
    Rcpp::IntegerVector output(*(context->S->ncol));
    int numVals = *(context->S->ncol);
    int i;
    for (i = 0; i < numVals; i++)
    {
        output[i] = (context->A0 -> R0)[i];
    }
    return(output);
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

Rcpp::NumericVector spatialSEIRInterface::getGenerationMatrix(int t)
{
    int nLoc = *(context->S->ncol);
    int numVals =nLoc*nLoc; 
    Rcpp::NumericMatrix output(nLoc, nLoc);
    double* input;
    int i;
    if (*(context -> isPopulated))
    {
        input = (context -> calculateG(t));
        for (i = 0; i < numVals; i++)
        {
            output[i] = input[i];
        }
        delete[] input;
        return(output);
    }
    Rcpp::Rcout << "Model context isn't populated\n";
    return(output);
}

Rcpp::NumericVector spatialSEIRInterface::getIntegratedGenerationMatrix(int t)
{
    int nLoc = *(context->S->ncol);
    int numVals =nLoc*nLoc; 
    Rcpp::NumericMatrix output(nLoc, nLoc);
    double* input;
    int i;
    if (*(context -> isPopulated))
    {
        input = (context -> calculateIntegratedG(t));
        for (i = 0; i < numVals; i++)
        {
            output[i] = input[i];
        }
        delete[] input;
        return(output);
    }
    Rcpp::Rcout << "Model context isn't populated\n";
    return(output);
}

Rcpp::NumericVector spatialSEIRInterface::getP_RS()
{
    Rcpp::NumericVector output(*(context->S->nrow));
    int i;
    int numVals = (*(context->S->nrow));
    for (i = 0; i < numVals; i++) 
    {
        output[i] = (context->p_rs)[i]; 
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
    int i;
    int numVals = (*(context -> S -> nrow));
    Rcpp::NumericVector output(numVals);
    for (i = 0; i < numVals; i ++) 
    {
        output[i] = (context->p_ei)[i]; 
    }
    return(output);
}
Rcpp::NumericVector spatialSEIRInterface::getP_IR()
{
    int i;
    int numVals = (*(context -> S -> nrow));
    Rcpp::NumericVector output(numVals);
    for (i = 0; i < numVals; i ++) 
    {
        output[i] = (context->p_ir)[i]; 
    }
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
                     SEXP Y_,
                     SEXP Sstar, 
                     SEXP Estar, 
                     SEXP Istar, 
                     SEXP Rstar, 
                     SEXP offset_,
                     SEXP X_,
                     SEXP Z_,
                     SEXP X_pRS_,
                     SEXP DistMat_,
                     SEXP rho_,
                     SEXP priorAlpha_pEI_,
                     SEXP priorBeta_pEI_,
                     SEXP priorAlpha_pIR_,
                     SEXP priorBeta_pIR_,
                     SEXP beta_,
                     SEXP betaPriorPrecision_,
                     SEXP betaPrs_,
                     SEXP betaPrsPriorPrecision_,
                     SEXP gamma_ei_,
                     SEXP gamma_ir_,
                     SEXP N_,
                     SEXP outFile,
                     SEXP iterationStride,
                     SEXP steadyStateConstraintPrecision_,
                     SEXP verboseFlag,
                     SEXP debugFlag,
                     SEXP sliceWidths,
                     SEXP reinfectionMode,
                     SEXP scaleDistanceMode_)
{
    int err = 0;
    //Deal with the data conversion from R to c++
    Rcpp::IntegerVector compartmentDimensions(compMatDim);
    Rcpp::IntegerVector covariateDimensions_x(xDim);
    Rcpp::IntegerVector covariateDimensions_z(zDim);
    Rcpp::IntegerVector covariateDimension_pRS_x(xPrsDim);
    Rcpp::IntegerVector S0(S0_);
    Rcpp::IntegerVector E0(E0_);
    Rcpp::IntegerVector I0(I0_);
    Rcpp::IntegerVector R0(R0_);

    Rcpp::IntegerVector Y(Y_);
    Rcpp::IntegerVector S_star(Sstar);
    Rcpp::IntegerVector E_star(Estar);
    Rcpp::IntegerVector I_star(Istar);
    Rcpp::IntegerVector R_star(Rstar);

    Rcpp::NumericVector offset(offset_);

    Rcpp::NumericVector X(X_);
    Rcpp::NumericVector Z(Z_);
    Rcpp::NumericVector X_pRS(X_pRS_);
    Rcpp::NumericVector DistMat(DistMat_);

    Rcpp::NumericVector rho(rho_);


    Rcpp::NumericVector priorAlpha_pEI(priorAlpha_pEI_);
    Rcpp::NumericVector priorBeta_pEI(priorBeta_pEI_);
    Rcpp::NumericVector priorAlpha_pIR(priorAlpha_pIR_);
    Rcpp::NumericVector priorBeta_pIR(priorBeta_pIR_);


    Rcpp::NumericVector beta(beta_);
    Rcpp::NumericVector betaPriorPrecision(betaPriorPrecision_);
    Rcpp::NumericVector betaPrs(betaPrs_);
    Rcpp::NumericVector betaPrsPriorPrecision(betaPrsPriorPrecision_);
    Rcpp::NumericVector gamma_ei(gamma_ei_);
    Rcpp::NumericVector gamma_ir(gamma_ir_);
    Rcpp::IntegerVector N(N_);

    Rcpp::NumericVector steadyStateConstraintPrecision(steadyStateConstraintPrecision_);
    Rcpp::NumericVector sliceParams(sliceWidths);

    Rcpp::IntegerVector reinfectMode(reinfectionMode);
    Rcpp::IntegerVector scaleDistanceMode(scaleDistanceMode_);

    chainOutputFile = new std::string(); 
    *chainOutputFile = Rcpp::as<std::string>(outFile);

    Rcpp::IntegerVector vFlag(verboseFlag);
    Rcpp::IntegerVector dFlag(debugFlag);
    *verbose = vFlag[0];
    *debug = dFlag[0];

    Rcpp::IntegerVector chainStride(iterationStride);
    

    try
    {

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
        if (S0.size() != compartmentDimensions[1])
        {
            Rcpp::Rcout << "Invalid S0 Compartment Size!\n";
            throw(-1);
        }
        if (E0.size() != compartmentDimensions[1])
        {
            Rcpp::Rcout << "Invalid E_star0 Compartment Size!\n";
            throw(-1);
        }
        if (I0.size() != compartmentDimensions[1])
        {
            Rcpp::Rcout << "Invalid I_star0 Compartment Size!\n";
            throw(-1);
        }
        if (R0.size() != compartmentDimensions[1])
        {
            Rcpp::Rcout << "Invalid R_star0 Compartment Size!\n";
            throw(-1);
        }

        if (Y.size() != compartmentSize)
        {
            Rcpp::Rcout << "Invalid Y Size!\n";
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
        if (offset.size() != compartmentDimensions[0])
        {
            Rcpp::Rcout << "Invalid Offset Size: " << offset.size() << "\n";
            throw(-1);
        }
        if (N.size() != compartmentSize)
        {
            Rcpp::Rcout << "Invalid N Compartment Size!\n";
            throw(-1);
        }
        if ((X_pRS.size() % compartmentDimensions[0]) != 0 && reinfectMode[0] <= 2)
        {
            Rcpp::Rcout << "Invalid X_pRS size.\n";
            Rcpp::Rcout << "Size: " << X_pRS.size() << ", Number of Time Points: " << compartmentDimensions[0] << "\n";
        }

        if (sliceParams.size() != 10)
        {
            Rcpp::Rcout << "Slice sampling parameters must be of length 10: S*,E*,R*,S0,I0,beta,betaPrs,rho,gamma_ei,gamma_ir\n";
            throw(-1);
        }
        if (reinfectMode[0] > 2)
        {
            int maxItr = (compartmentDimensions[0]*compartmentDimensions[1]); 
            int i;
            for (i = 0; i < maxItr; i++)
            {
                if (S_star[i] != 0)
                {
                    Rcpp::Rcout << "Error: reinfectionMode indicates that no reinfection shoul occur, but nonzero S_star provided\n";
                    throw(-1);
                }
            }
        }
        if (scaleDistanceMode[0] != 0 && scaleDistanceMode[0] != 1)
        {
            Rcpp::Rcout << "scaleDistanceMode should be true or false until more options are implemented.\n";
            Rcpp::Rcout << "   mode given: " << scaleDistanceMode[0] << "\n";
            throw(-1);
        }
    }
    catch(int e)
    {
        Rcpp::Rcout << "Errors Encountered, exiting.\n";
        delete chainOutputFile;
        return -1;
    }





    Rcpp::Rcout << "Building Model.\n   Number of Locations: " << compartmentDimensions[1] 
        << "\n";
    Rcpp::Rcout << "   Number of Time Points: " << compartmentDimensions[0] 
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
    sliceParameters sliceParamStruct;
    modelConfiguration modelConfig;
    modelConfig.reinfectionMode = reinfectMode[0];
    modelConfig.compartmentSamplingMode = COMPARTMENT_METROPOLIS_SAMPLER;
    modelConfig.parameterSamplingMode = PARAMETER_JOINT_METROPOLIS_SAMPLER;
    modelConfig.indexLength = std::floor(0.25*compartmentDimensions[0]*compartmentDimensions[1]); // Update 25% per iteration. 
    modelConfig.useDecorrelation = 0;
    Rcpp::Rcout << "Setting index length to be: " << (modelConfig.indexLength) << "\n";

    sliceParamStruct.S_starWidth = &sliceParams[0];
    sliceParamStruct.E_starWidth = &sliceParams[1];
    sliceParamStruct.R_starWidth = &sliceParams[2];
    sliceParamStruct.S0Width = &sliceParams[3];
    sliceParamStruct.I0Width = &sliceParams[4];
    sliceParamStruct.betaWidth = &sliceParams[5];
    sliceParamStruct.betaPrsWidth = &sliceParams[6];
    sliceParamStruct.rhoWidth = &sliceParams[7];
    sliceParamStruct.gammaEiWidth = &sliceParams[8];
    sliceParamStruct.gammaIrWidth = &sliceParams[9];

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


    priorControl priorValues;
    priorValues.betaPriorPrecision = betaPriorPrecision[0];
    priorValues.P_EI_priorAlpha = priorAlpha_pEI[0];
    priorValues.P_EI_priorBeta = priorBeta_pEI[0];
    priorValues.P_IR_priorAlpha = priorAlpha_pIR[0];
    priorValues.P_IR_priorBeta = priorBeta_pIR[0];
    priorValues.betaPrsPriorPrecision = betaPrsPriorPrecision[0];

    // Gather information for the creation of the distance matrices

    double phi = 60*60*4.0;
    distanceArgs rawDistArgs; scaledDistanceArgs scaledDistArgs;
    rawDistArgs.inData = DistMat.begin(); 
    rawDistArgs.dim = &compartmentDimensions[1];
    scaledDistArgs.phi = &phi; 
    scaledDistArgs.inData = DistMat.begin();
    scaledDistArgs.dim = &compartmentDimensions[1];
    scaledDistArgs.mode = scaleDistanceMode[0];

    // Create the InitData object 
    InitData A0;
    A0.populate(S0.begin(),E0.begin(),I0.begin(),R0.begin()
            ,&compartmentDimensions[1]);

    //Rcpp::Rcout << compartmentDimensions[0] << " " << compartmentDimensions[1] << "\n";
    //Rcpp::Rcout << (xArgs.inData_x)[1] << "\n";
    context -> populate(&A0, &xArgs, &xPrsArgs, offset.begin(), Y.begin(), &S_starArgs, &E_starArgs, &I_starArgs, 
                        &R_starArgs, &rawDistArgs,&scaledDistArgs,
                        rho.begin(),beta.begin(),gamma_ei.begin(), gamma_ir.begin(),
                        betaPrs.begin(),N.begin(),&sliceParamStruct, &priorValues,
                        modelConfig);

    // Set up output stream
    context -> fileProvider -> populate(context, chainOutputFile,
            (int*) chainStride.begin());

    return(err);
}


spatialSEIRInterface::~spatialSEIRInterface()
{   
    // Context handles the complicated cleanup
    delete verbose;
    delete debug;
    delete context;
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
    .method("setPredictionTraces", &spatialSEIRInterface::setPredictionTraces)
    .method("setTrace", &spatialSEIRInterface::setTrace)
    .method("setTrace", &spatialSEIRInterface::setTrace2)
    .method("setDevice", &spatialSEIRInterface::setDevice)
    .method("calculateS", &spatialSEIRInterface::calculateS)
    .method("calculateE", &spatialSEIRInterface::calculateE)
    .method("calculateI", &spatialSEIRInterface::calculateI)
    .method("calculateR", &spatialSEIRInterface::calculateR)
    .method("calculateP_RS", &spatialSEIRInterface::calculateP_RS)
    .method("calculateP_SE", &spatialSEIRInterface::calculateP_SE)
    .method("calculateP_SE", &spatialSEIRInterface::calculateP_SE2)
    .method("calculateP_SE_OCL", &spatialSEIRInterface::calculateP_SE_OCL)
    .method("estimateR0", &spatialSEIRInterface::estimateR03)
    .method("estimateR0", &spatialSEIRInterface::estimateR02)
    .method("estimateR0", &spatialSEIRInterface::estimateR0)
    .method("printAcceptanceRates", &spatialSEIRInterface::printAcceptanceRates)
    .method("printOCLSummary", &spatialSEIRInterface::printOCLSummary)
    .method("printSamplingParameters", &spatialSEIRInterface::printSamplingParameters)
    .method("updateSamplingParameters", &spatialSEIRInterface::updateSamplingParameters)
    .method("getGenerationMatrix", &spatialSEIRInterface::getGenerationMatrix)
    .method("getIntegratedGenerationMatrix", &spatialSEIRInterface::getIntegratedGenerationMatrix)
    .method("standardizeDistanceMatrix", &spatialSEIRInterface::standardizeDistanceMatrix)
    .property("S", &spatialSEIRInterface::getS, "Susceptible Compartment Matrix")
    .property("E", &spatialSEIRInterface::getE, "Exposed Compartment Matrix")
    .property("I", &spatialSEIRInterface::getI, "Infectious Compartment Matrix")
    .property("R", &spatialSEIRInterface::getR, "Removed Compartment Matrix")
    .property("S0", &spatialSEIRInterface::getS0, "Initial Susceptible Compartment Matrix")
    .property("E0", &spatialSEIRInterface::getE0, "Initial Exposed Compartment Matrix")
    .property("I0", &spatialSEIRInterface::getI0, "Initial Infectious Compartment Matrix")
    .property("R0", &spatialSEIRInterface::getR0, "Initial Removed Compartment Matrix")
    .property("S_star", &spatialSEIRInterface::getS_star, "Removed to Susceptible Transition Matrix")
    .property("E_star", &spatialSEIRInterface::getE_star, "Susceptible to Exposed Transition Matrix")
    .property("I_star", &spatialSEIRInterface::getI_star, "Exposed to Infectious Transition Matrix")
    .property("R_star", &spatialSEIRInterface::getR_star, "Infectious to Removed Transition Matrix")
    .property("p_se", &spatialSEIRInterface::getP_SE, "Exposure Probability Matrix")
    .property("p_ei", &spatialSEIRInterface::getP_EI, "E to I Transition Probability")
    .property("p_ir", &spatialSEIRInterface::getP_IR, "I to R Transition Probability")
    .property("p_rs", &spatialSEIRInterface::getP_RS, "R-S Transition Probability Vector")
    .property("beta", &spatialSEIRInterface::getBeta, "Exposure Process Regression Parameters")
    .property("betaP_RS", &spatialSEIRInterface::getBetaP_RS, "R-S Transition Process Regression Parameters")
    .property("rho", &spatialSEIRInterface::getRho, "Spatial Dependence Term")
    .property("parameterSamplingMode", &spatialSEIRInterface::getParameterSamplingMode, 
            &spatialSEIRInterface::setParameterSamplingMode, "Type of sampler used for non-compartment parameters.")
    .property("compartmentSamplingMode", &spatialSEIRInterface::getCompartmentSamplingMode, 
            &spatialSEIRInterface::setCompartmentSamplingMode, "Type of sampler used for disease compartments.")
    .property("debug", &spatialSEIRInterface::getDebug, &spatialSEIRInterface::setDebug, "Show debug level output?")
    .property("verbose", &spatialSEIRInterface::getVerbose, &spatialSEIRInterface::setVerbose, "Show verbose level output?")
    .property("useDecorrelation", &spatialSEIRInterface::getUseDecorrelation, &spatialSEIRInterface::setUseDecorrelation, "Use decorrelation sampling?")
    ;

}

