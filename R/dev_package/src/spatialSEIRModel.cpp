#include <Rcpp.h>
#include <Rmath.h>
#include <spatialSEIRModel.hpp>
#include <cmath>
#include <dataModel.hpp>
#include <exposureModel.hpp>
#include <reinfectionModel.hpp>
#include <distanceModel.hpp>
#include <transitionPriors.hpp>
#include <initialValueContainer.hpp>
#include <samplingControl.hpp>
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

int spatialSEIRModel::getUseDecorrelation()
{
    if (*(context -> isPopulated))
    {
        return((context -> config -> useDecorrelation));
    }
    Rcpp::Rcout << "Context Not populated\n";
    return(-1);
}

void spatialSEIRModel::setUseDecorrelation(int val)
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

int spatialSEIRModel::setRandomSeed(int seedVal)
{
    if (*(context -> isPopulated))
    {
        context -> setRandomSeed(seedVal);
        return(0);
    }
    Rcpp::Rcout << "Context Not populated\n";
    return(-1);
}
int spatialSEIRModel::simulate(int iters)
{
    if (*(context -> isPopulated))
    {
        context -> runSimulation(iters,*(verbose),*(debug));
        return(0);
    }
    Rcpp::Rcout << "Context Not populated\n";

}

void spatialSEIRModel::setDevice(int platformId, int deviceId)
{
    if (*(context -> isPopulated))
    {
        (context -> oclProvider -> setDevice(platformId, deviceId));
        return;
    }
    Rcpp::Rcout << "ModelContext is not populated.\n";
    return;
}

void spatialSEIRModel::setCompartmentSamplingMode(int mode)
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
        else if (mode == COMPARTMENT_BINOM_IDX_METROPOLIS_SAMPLER)
        {
            Rcpp::Rcout << "Setting compartment sampling mode to indexed Metropolis with chain binomial proposal.\n";
        }
        else
        {
            Rcpp::Rcout << "Error: mode must be one of the following: \n"
                        << COMPARTMENT_METROPOLIS_SAMPLER << ": Metropolis \n   "
                        << COMPARTMENT_METROPOLIS_SAMPLER_OCL << ": Metropolis with OpenCL\n   "
                        << COMPARTMENT_IDX_METROPOLIS_SAMPLER << ": Indexed Metropolis\n   "
                        << COMPARTMENT_BINOM_PROPOSAL_METROPOLIS_SAMPLER << ": Metropolis w/ chain binomial proposal.\n   " 
                        << COMPARTMENT_BINOM_IDX_METROPOLIS_SAMPLER << ": Indexed Metropolis w/ chain binomial proposal.\n   " 
                        << COMPARTMENT_BINOM_PROPOSAL_SLICE_SAMPLER << ": Slice w/ chain binomial proposal.\n   " 
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

int spatialSEIRModel::getCompartmentSamplingMode()
{
    if (*(context -> isPopulated))
    {
        return((context -> getCompartmentSamplingMode()));
    }
    return(-1);
}

void spatialSEIRModel::setParameterSamplingMode(int mode)
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

int spatialSEIRModel::getParameterSamplingMode()
{
    if (*(context -> isPopulated))
    {
        return((context -> getParameterSamplingMode()));
    }
    return(-1);
}



void spatialSEIRModel::standardizeDistanceMatrices()
{
    if (*(context -> isPopulated))
    {
       if (*(context -> numIterations) != 0)
       {
           Rcpp::Rcout << "Can't change distance matrix once sampling has begun.\n";
           return;
       }
       unsigned int k;
       for (k = 0; k < context -> scaledDistMatrices -> size(); k++)
        {
           (*(context -> scaledDistMatrices))[k] -> makeRowStochastic(); 
        }
    }
    else
    {
        Rcpp::Rcout << "No distance matrix to standardize.\n";
    }
}


int spatialSEIRModel::setTrace(int locationIndex)
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
int spatialSEIRModel::setTrace2(int locationIndex, int timeIndex)
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

void spatialSEIRModel::setPredictionTraces()
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

int spatialSEIRModel::printDebugInfo()
{
    if (*(context -> isPopulated))
    {
        Rcpp::Rcout << "S Dimensions: ";
        Rcpp::Rcout << *(context -> S -> nrow) << ", " << *(context -> S -> ncol) << "\n";
        return(0);
    }
    Rcpp::Rcout << "Context not populated\n";
    return(-1);
}

int spatialSEIRModel::calculateS()
{
    if (*(context -> isPopulated))
    {
        (context -> calculateS_CPU());
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate S on non-populated ModelContext.\n";
    return(-1);
}
int spatialSEIRModel::calculateE()
{
    if (*(context -> isPopulated))
    {
        (context -> calculateE_CPU());
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate E on non-populated ModelContext.\n";
    return(-1);
}
int spatialSEIRModel::calculateI()
{
    if (*(context -> isPopulated))
    {
        (context -> calculateI_CPU());
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate I on non-populated ModelContext.\n";
    return(-1);
}
int spatialSEIRModel::calculateR()
{
    if (*(context -> isPopulated))
    {
        (context -> calculateR_CPU());
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate R on non-populated ModelContext.\n";
    return(-1);
}
int spatialSEIRModel::calculateP_SE() 
{
    if (*(context -> isPopulated))
    {
        (context -> calculateP_SE_CPU());
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate P_SE on non-populated ModelContext.\n";
    return(-1);
}
int spatialSEIRModel::calculateP_SE2(int i, int j) 
{
    if (*(context -> isPopulated))
    {
        (context -> calculateP_SE_CPU(i,j));
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate P_SE on non-populated ModelContext.\n";
    return(-1);
}
int spatialSEIRModel::calculateP_SE_OCL() 
{
    if (*(context -> isPopulated))
    {
        (context -> calculateP_SE_OCL());
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate P_SE on non-populated ModelContext.\n";
    return(-1);
}

double spatialSEIRModel::estimateR0()
{
    try
    {
        if (*(context -> isPopulated))
        {
            return((context -> estimateR0()));
        }
        Rcpp::Rcout << "Attempt to estimate R0 on a non-populated ModelContext.\n";
        return(-1.0);
    }
    catch (int e)
    {
        Rcpp::Rcout << "Error: " << e << "\n";
        return(-1.0);
    }
}

double spatialSEIRModel::estimateEffectiveR0()
{
    try
    {
        if (*(context -> isPopulated))
        {
            return((context -> estimateEffectiveR0()));
        }
        Rcpp::Rcout << "Attempt to estimate R0 on a non-populated ModelContext.\n";
        return(-1.0);
    }
    catch (int e)
    {
        Rcpp::Rcout << "Error: " << e << "\n";
        return(-1.0);
    }
}

Rcpp::NumericVector spatialSEIRModel::estimateR02(int t)
{
    try
    {
        if (*(context -> isPopulated))
        {
            int nLoc = *(context -> S -> ncol);
            double* inVec = (context -> estimateR0(t));
            Rcpp::NumericVector outVec (nLoc);
            memcpy(outVec.begin(), inVec, nLoc*sizeof(double));
            return(outVec);
        }
        Rcpp::Rcout << "Attempt to estimate R0 on a non-populated ModelContext.\n";
        return(-1.0);
    }
    catch (int e)
    {
        Rcpp::Rcout << "Error: " << e << "\n";
        return(-1.0);
    }

}

Rcpp::NumericVector spatialSEIRModel::estimateEffectiveR02(int t)
{
    try
    {
        if (*(context -> isPopulated))
        {
            int nLoc = *(context -> S -> ncol);
            double* inVec = (context -> estimateEffectiveR0(t));
            Rcpp::NumericVector outVec (nLoc);
            memcpy(outVec.begin(), inVec, nLoc*sizeof(double));
            return(outVec);
        }
        Rcpp::Rcout << "Attempt to estimate R0 on a non-populated ModelContext.\n";
        return(-1.0);
    }
    catch (int e)
    {
        Rcpp::Rcout << "Error: " << e << "\n";
        return(-1.0);
    }

}

int spatialSEIRModel::calculateP_RS()
{
    if (*(context -> isPopulated))
    {
        (context -> calculateP_RS_CPU());
        return(0);
    }
    Rcpp::Rcout << "Attept to calculate P_RS on non-populated ModelContext.\n";
    return(-1);
}

void spatialSEIRModel::updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange)
{
    if ((*(context -> numIterations)) == 0 ||  !(*(context -> isPopulated)))
    {
        Rcpp::Rcout << "No samples drawn.\n";
        return;
    }
    (context -> updateSamplingParameters(desiredRatio, targetWidth, proportionChange));
}

void spatialSEIRModel::printOCLSummary()
{
    if (*(context -> isPopulated))
    {
        (context -> oclProvider -> printSummary());
        return;
    }
    Rcpp::Rcout << "ModelContext has not been populated.\n";
}

void spatialSEIRModel::printSamplingParameters()
{
    if (*(context -> isPopulated))
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
        return;
    }
    Rcpp::Rcout << "Context not populated\n";
}

void spatialSEIRModel::printAcceptanceRates()
{
    if (*(context -> isPopulated))
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
            Rcpp::Rcout << "rho:      ";
            for (i = 0; i < *(context -> rho_fc -> varLen); i++)
            {

                Rcpp::Rcout << ((context -> rho_fc -> accepted)[i]*1.0)/
                                                  (*(context -> rho_fc -> samples)) 
                                          << ", "; 

            }
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
        return;
    }
    Rcpp::Rcout << "Context not populated\n";
}

Rcpp::IntegerVector spatialSEIRModel::getS0()
{
    if (*(context -> isPopulated))
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
    Rcpp::Rcout << "Context not populated\n";
    Rcpp::IntegerVector output(1);
    output[0] = -1;
    return(output);
}
Rcpp::IntegerVector spatialSEIRModel::getE0()
{
    if (*(context -> isPopulated))
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
    Rcpp::Rcout << "Context not populated\n";
    Rcpp::IntegerVector output(1);
    output[0] = -1;
    return(output);

}
Rcpp::IntegerVector spatialSEIRModel::getI0()
{
    if (*(context -> isPopulated))
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
    Rcpp::Rcout << "Context not populated\n";
    Rcpp::IntegerVector output(1);
    output[0] = -1;
    return(output);
}
Rcpp::IntegerVector spatialSEIRModel::getR0()
{
    if (*(context -> isPopulated))
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
    Rcpp::Rcout << "Context not populated\n";
    Rcpp::IntegerVector output(1);
    output[0] = -1;
    return(output);

}

Rcpp::IntegerMatrix spatialSEIRModel::getS()
{
    if (*(context -> isPopulated))
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
    Rcpp::Rcout << "Context not populated\n";
    Rcpp::IntegerMatrix output(1,1);
    output[0] = -1;
    return(output);
}
Rcpp::IntegerMatrix spatialSEIRModel::getE()
{
    if (*(context -> isPopulated))
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
    Rcpp::Rcout << "Context not populated\n";
    Rcpp::IntegerMatrix output(1,1);
    output[0] = -1;
    return(output);
}
Rcpp::IntegerMatrix spatialSEIRModel::getI()
{
    if (*(context -> isPopulated))
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
    Rcpp::Rcout << "Context not populated\n";
    Rcpp::IntegerMatrix output(1,1);
    output[0] = -1;
    return(output);

}
Rcpp::IntegerMatrix spatialSEIRModel::getR()
{
    if (*(context -> isPopulated))
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
    Rcpp::Rcout << "Context not populated\n";
    Rcpp::IntegerMatrix output(1,1);
    output[0] = -1;
    return(output);
}

Rcpp::IntegerMatrix spatialSEIRModel::getY()
{
    if (*(context -> isPopulated))
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
    Rcpp::Rcout << "Context Not Populated\n";
    Rcpp::IntegerMatrix err(0,0);
    err[0] = -1;

}

Rcpp::IntegerMatrix spatialSEIRModel::getS_star()
{
    if (*(context -> isPopulated))
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
    Rcpp::Rcout << "Context not populated\n";
    Rcpp::IntegerMatrix output(1,1);
    output[0] = -1;
    return(output);
}

Rcpp::IntegerMatrix spatialSEIRModel::getE_star()
{
    if (*(context -> isPopulated))
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
    Rcpp::Rcout << "Context not populated\n";
    Rcpp::IntegerMatrix output(1,1);
    output[0] = -1;
    return(output);
}

Rcpp::IntegerMatrix spatialSEIRModel::getI_star()
{
    if (*(context -> isPopulated))
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
    Rcpp::Rcout << "Context not populated\n";
    Rcpp::IntegerMatrix output(1,1);
    output[0] = -1;
    return(output);
}

Rcpp::IntegerMatrix spatialSEIRModel::getR_star()
{
    if (*(context -> isPopulated))
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
    Rcpp::Rcout << "Context not populated\n";
    Rcpp::IntegerMatrix output(1,1);
    output[0] = -1;
    return(output);
}

Rcpp::NumericMatrix spatialSEIRModel::getP_SE()
{
    if (*(context -> isPopulated))
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
    Rcpp::Rcout << "Context not populated\n";
    Rcpp::NumericMatrix output(1,1);
    output[0] = -1;
    return(output);
}

Rcpp::NumericVector spatialSEIRModel::getGenerationMatrix(int t)
{
    if (!*(context -> isPopulated))
    {
        Rcpp::Rcout << "Context not populated\n";
        throw(-1);
    }
    try
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
    catch (int e)
    {
        Rcpp::Rcout << "Error: " << e << "\n";
        Rcpp::NumericMatrix out(1,1);
        out[0] = -1;
        return(out);
    }
}

Rcpp::NumericVector spatialSEIRModel::getIntegratedGenerationMatrix(int t)
{
    if (!*(context -> isPopulated))
    {
        Rcpp::Rcout << "Context not populated\n";
        throw(-1);
    }
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

Rcpp::NumericVector spatialSEIRModel::getP_RS()
{
    if (!*(context -> isPopulated))
    {
        Rcpp::Rcout << "Context not populated\n";
        throw(-1);
    }

    Rcpp::NumericVector output(*(context->S->nrow));
    int i;
    int numVals = (*(context->S->nrow));
    for (i = 0; i < numVals; i++) 
    {
        output[i] = (context->p_rs)[i]; 
    }
    return(output);
}

Rcpp::NumericVector spatialSEIRModel::getBeta()
{
    if (!*(context -> isPopulated))
    {
        Rcpp::Rcout << "Context not populated\n";
        throw(-1);
    }

    Rcpp::NumericVector output(((*(context->X->ncol_x)) + (*(context->X->ncol_z))));
    int i;
    int numVals = (((*(context->X->ncol_x)) + (*(context->X->ncol_z))));
    for (i = 0; i < numVals; i++) 
    {
        output[i] = (context->beta)[i]; 
    }
    return(output);

}
Rcpp::NumericVector spatialSEIRModel::getBetaP_RS()
{
    if (!*(context -> isPopulated))
    {
        Rcpp::Rcout << "Context not populated\n";
        throw(-1);
    }

    Rcpp::NumericVector output(*(context->X_pRS->ncol_x));
    int i;
    int numVals = (*(context->X_pRS->ncol_x));
    for (i = 0; i < numVals; i++) 
    {
        output[i] = (context->betaPrs)[i]; 
    }
    return(output);
}
Rcpp::NumericVector spatialSEIRModel::getP_EI()
{
    if (!*(context -> isPopulated))
    {
        Rcpp::Rcout << "Context not populated\n";
        throw(-1);
    }

    int i;
    int numVals = (*(context -> S -> nrow));
    Rcpp::NumericVector output(numVals);
    for (i = 0; i < numVals; i ++) 
    {
        output[i] = (context->p_ei)[i]; 
    }
    return(output);
}
Rcpp::NumericVector spatialSEIRModel::getP_IR()
{
    if (!*(context -> isPopulated))
    {
        Rcpp::Rcout << "Context not populated\n";
        throw(-1);
    }

    int i;
    int numVals = (*(context -> S -> nrow));
    Rcpp::NumericVector output(numVals);
    for (i = 0; i < numVals; i ++) 
    {
        output[i] = (context->p_ir)[i]; 
    }
    return(output);
}

Rcpp::NumericVector spatialSEIRModel::getRho()
{
    if (*(context -> isPopulated))
    {
        int i;
        int nRho = (context -> scaledDistMatrices -> size()); 
        Rcpp::Rcout << nRho << "\n";
        Rcpp::NumericVector output(nRho);
        for (i = 0; i < nRho; i ++)
        {
           output[i] = (context -> rho)[i]; 
        }
        return(output);
    }
    Rcpp::Rcout << "Context Not Populated\n";
    Rcpp::NumericVector output(1);
    output[0] = -1.0;
}

Rcpp::NumericVector spatialSEIRModel::getPhi()
{
    Rcpp::NumericVector output(1);
    if (*(context -> isPopulated))
    {
        output[0] = *(context->phi); 
        return(output);
    }
    Rcpp::Rcout << "Context Not Populated\n";
    output[0] = -1.0;
}



int spatialSEIRModel::getDebug()
{
    int out;
    out = *debug;
    return(out);    
}
void spatialSEIRModel::setDebug(int debug_)
{
   *debug = debug_;  
}
int spatialSEIRModel::getVerbose()
{
    int output;
    output= *verbose;
    return(output);
}
void spatialSEIRModel::setVerbose(int verbose_)
{
   *verbose = verbose_; 
}

spatialSEIRModel::spatialSEIRModel(SEXP outFileName)
{
    // Create the empty ModelContext object  
    context = new ModelContext();
    verbose = new int();
    debug = new int(); 
    *verbose = 0;
    *debug = 0;
    chainOutputFile = new std::string(); 
    *chainOutputFile = Rcpp::as<std::string>(outFileName);
}

int spatialSEIRModel::buildSpatialSEIRModel(dataModel& dataModel_,
                          exposureModel& exposureModel_,
                          reinfectionModel& reinfectionModel_,
                          distanceModel& distanceModel_,
                          transitionPriors& transitionPriors_,
                          initialValueContainer& initialValueContainer_,
                          samplingControl& samplingControl_)
{
    int err = 0;
    
    dataModel* dataModelInstance = &dataModel_;
    exposureModel* exposureModelInstance = &exposureModel_;
    reinfectionModel* reinfectionModelInstance = &reinfectionModel_;
    distanceModel* distanceModelInstance = &distanceModel_;
    transitionPriors* transitionPriorsInstance = &transitionPriors_;
    initialValueContainer* initialValueContainerInstance = &initialValueContainer_;
    samplingControl* samplingControlInstance = &samplingControl_;


    if (*(dataModelInstance -> nLoc) != (exposureModelInstance -> xDim)[0])
    {
        Rcpp::Rcout << "Exposure model and data model imply different number of locations\n";
        Rcpp::Rcout << "Exposure model: " << (exposureModelInstance -> xDim)[0] << "\n";
        Rcpp::Rcout << "Data model: " << *(dataModelInstance -> nLoc) << "\n";
        throw(-1);
    }
    if (*(dataModelInstance -> nTpt) != ((exposureModelInstance -> zDim)[0])/((exposureModelInstance -> xDim)[0]))
    {
        Rcpp::Rcout << "Exposure model and data model imply different number of time points\n";
        throw(-1);
    }
    if (*(dataModelInstance -> nLoc) != (*(distanceModelInstance -> numLocations)))
    {
        Rcpp::Rcout << "Data model and distance model imply different number of locations\n";
        throw(-1);
    }
    if ((*(dataModelInstance -> nLoc) != (initialValueContainerInstance -> compMatDim)[1]) || 
        (*(dataModelInstance -> nTpt) != (initialValueContainerInstance -> compMatDim)[0]))
    {
        Rcpp::Rcout << "Data model and initial value container have different dimensions\n";
        Rcpp::Rcout << "Data Model: " << *(dataModelInstance -> nTpt) << ", " 
            << *(dataModelInstance -> nLoc) << "\n";
        Rcpp::Rcout << "Initial Vals: " << (initialValueContainerInstance -> compMatDim)[0] << ", "
            << (initialValueContainerInstance -> compMatDim)[1] << "\n";
        throw(-1);
    }
    if (*(reinfectionModelInstance -> reinfectionMode) == 3)
    {
        // No reinfection
    }
    else
    {
        if ((reinfectionModelInstance -> xDim)[0] != *(dataModelInstance -> nTpt))
        {
            Rcpp::Rcout << "Reinfection and data mode time points differ.\n";
            throw(-1);
        }
    }
    if (*(transitionPriorsInstance -> gamma_ei) < 0 || 
        *(transitionPriorsInstance -> gamma_ir) < 0)
    {
        Rcpp::Rcout << "Transition priors haven't been populated (or have invalid values)\n";
        throw(-1);
    }
    int i;
    if (*(reinfectionModelInstance -> reinfectionMode) > 2)
    {
        int maxItr = ((*(dataModelInstance -> compartmentDimensions))[0]
                     *(*(dataModelInstance -> compartmentDimensions))[1]);
        for (i = 0; i < maxItr; i++)
        {
            if ((initialValueContainerInstance -> S_star)[i] != 0)
            {
                Rcpp::Rcout << "Error: reinfectionMode indicates that no reinfection should occur, but nonzero S_star provided\n";
                throw(-1);
            }
        }
        // Create appropriate dummy data. 
        reinfectionModelInstance -> buildDummyReinfectionModel((*(dataModelInstance -> compartmentDimensions))[0]);
    }

    if (*(dataModelInstance -> dataModelType) != 0 && *(dataModelInstance ->setMode) < 0)
    {
        Rcpp::Rcout << "Non-identity data model requested, but prior parameters were not supplied.\n"; 
        throw(-1);
    }

    int* nTpt = &((*(dataModelInstance -> compartmentDimensions))[0]); 
    int* nLoc = &((*(dataModelInstance -> compartmentDimensions))[1]);
    
    Rcpp::Rcout << "Building Model.\n   Number of Locations: " << *nLoc
        << "\n";
    Rcpp::Rcout << "   Number of Time Points: " << *nTpt
        << "\n";

    
    int numDistMatrices = (distanceModelInstance -> getNumDistanceMatrices());
    Rcpp::NumericVector rho(numDistMatrices);
    double rhoSum = 0.0;
    for (i = 0; i < numDistMatrices; i++)
    {
        rho[i] = R::rgamma(0.5,0.5);
        rhoSum += rho[i];
    }
    if (rhoSum >= 0.5)
    {
        for (i = 0; i < numDistMatrices; i++)
        {
            rho[i]*=(1/rhoSum)*0.1;
        }
    }

    Rcpp::NumericVector phi(1);
    phi[0] = (dataModelInstance -> initialParameterValues)[0];

    double* sliceParams = (samplingControlInstance -> sliceWidths); 

    *verbose = *(samplingControlInstance -> verbose);
    *debug = *(samplingControlInstance -> debug);
    

    // Gather information for the creation of the 
    // covariate matrix
    covariateArgs xArgs;
    xArgs.inData_x = (exposureModelInstance -> X);
    xArgs.inData_z = (exposureModelInstance -> Z);
    xArgs.inRow_x = &((exposureModelInstance -> xDim)[0]);
    xArgs.inCol_x = &((exposureModelInstance -> xDim)[1]);
    xArgs.inRow_z = &((exposureModelInstance -> zDim)[0]);
    xArgs.inCol_z = &((exposureModelInstance -> zDim)[1]);

    covariateArgs xPrsArgs; 
    xPrsArgs.inData_x = (reinfectionModelInstance -> X);
    xPrsArgs.inData_z = NULL;
    xPrsArgs.inRow_x = &((reinfectionModelInstance -> xDim)[0]);
    xPrsArgs.inCol_x = &((reinfectionModelInstance -> xDim)[1]);
    // Clean this up, pass values instead. 
    int zeroVal = 0;
    xPrsArgs.inRow_z = &zeroVal;
    xPrsArgs.inCol_z = &zeroVal;


    // Gather information for the creation of the compartment matrices 
    
    compartmentArgs S_starArgs, E_starArgs, I_starArgs, R_starArgs;
    sliceParameters sliceParamStruct;
    modelConfiguration modelConfig;
    modelConfig.reinfectionMode = *(reinfectionModelInstance -> reinfectionMode);
    modelConfig.compartmentSamplingMode = COMPARTMENT_METROPOLIS_SAMPLER;
    modelConfig.parameterSamplingMode = PARAMETER_JOINT_METROPOLIS_SAMPLER;
    modelConfig.indexLength = std::floor(0.25*(*nTpt)*(*nLoc)); // Update 25% per iteration. 
    modelConfig.useDecorrelation = 0;
    modelConfig.dataModel = *(dataModelInstance -> dataModelType) ;

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
    sliceParamStruct.phiWidth = &sliceParams[10];

    S_starArgs.inData = (initialValueContainerInstance -> S_star);
    S_starArgs.inRow = nTpt;
    S_starArgs.inCol = nLoc;
    S_starArgs.steadyStateConstraintPrecision = (samplingControlInstance -> steadyStateConstraintPrecision)[0];

    E_starArgs.inData = (initialValueContainerInstance -> E_star);
    E_starArgs.inRow = nTpt;
    E_starArgs.inCol = nLoc;
    E_starArgs.steadyStateConstraintPrecision = (samplingControlInstance -> steadyStateConstraintPrecision)[0];

    I_starArgs.inData = (initialValueContainerInstance -> I_star);
    I_starArgs.inRow = nTpt;
    I_starArgs.inCol = nLoc;

    R_starArgs.inData = (initialValueContainerInstance -> R_star);
    R_starArgs.inRow = nTpt;
    R_starArgs.inCol = nLoc;
    R_starArgs.steadyStateConstraintPrecision = (samplingControlInstance -> steadyStateConstraintPrecision)[0];


    priorControl priorValues;
    priorValues.betaPriorPrecision = *(exposureModelInstance -> betaPriorPrecision);
    priorValues.P_EI_priorAlpha = (transitionPriorsInstance -> gamma_ei_params)[0];
    priorValues.P_EI_priorBeta = (transitionPriorsInstance -> gamma_ei_params)[1];
    priorValues.P_IR_priorAlpha = (transitionPriorsInstance -> gamma_ir_params)[0];
    priorValues.P_IR_priorBeta = (transitionPriorsInstance -> gamma_ir_params)[1];
    priorValues.betaPrsPriorPrecision = *(reinfectionModelInstance -> betaPriorPrecision);
    priorValues.Phi_priorAlpha = (dataModelInstance -> priorParameters)[0];
    priorValues.Phi_priorBeta = (dataModelInstance -> priorParameters)[1];
    

    // Create the InitData object 
    InitData A0;
    A0.populate((initialValueContainerInstance -> S0),
                (initialValueContainerInstance -> E0),
                (initialValueContainerInstance -> I0),
                (initialValueContainerInstance -> R0),
                nLoc);

    context -> populate(&A0, &xArgs, &xPrsArgs, (exposureModelInstance -> offset), (dataModelInstance -> Y), &S_starArgs, &E_starArgs, &I_starArgs, 
                        &R_starArgs, distanceModelInstance -> scaledDistArgs,
                        rho.begin(),phi.begin(),(exposureModelInstance -> beta),(transitionPriorsInstance -> gamma_ei), (transitionPriorsInstance -> gamma_ir),
                        (reinfectionModelInstance -> beta), (initialValueContainerInstance -> N),&sliceParamStruct, &priorValues,
                        modelConfig);

    // Set up output stream
    context -> fileProvider -> populate(context, chainOutputFile,(samplingControlInstance -> iterationStride));
    return(err);
}


spatialSEIRModel::~spatialSEIRModel()
{   
    // Context handles the complicated cleanup
    delete verbose;
    delete debug;
    delete context;
    delete chainOutputFile;
}


RCPP_MODULE(mod_spatialSEIRModel)
{
    using namespace Rcpp;
    class_<spatialSEIRModel>( "spatialSEIRModel" )

    .constructor<SEXP>()

    .method("buildSpatialSEIRModel", &spatialSEIRModel::buildSpatialSEIRModel)
    .method("printDebugInfo", &spatialSEIRModel::printDebugInfo)
    .method("setRandomSeed", &spatialSEIRModel::setRandomSeed)
    .method("simulate", &spatialSEIRModel::simulate)
    .method("setPredictionTraces", &spatialSEIRModel::setPredictionTraces)
    .method("setTrace", &spatialSEIRModel::setTrace)
    .method("setTrace", &spatialSEIRModel::setTrace2)
    .method("setDevice", &spatialSEIRModel::setDevice)
    .method("calculateS", &spatialSEIRModel::calculateS)
    .method("calculateE", &spatialSEIRModel::calculateE)
    .method("calculateI", &spatialSEIRModel::calculateI)
    .method("calculateR", &spatialSEIRModel::calculateR)
    .method("calculateP_RS", &spatialSEIRModel::calculateP_RS)
    .method("calculateP_SE", &spatialSEIRModel::calculateP_SE)
    .method("calculateP_SE", &spatialSEIRModel::calculateP_SE2)
    .method("calculateP_SE_OCL", &spatialSEIRModel::calculateP_SE_OCL)
    .method("estimateR0", &spatialSEIRModel::estimateR02)
    .method("estimateR0", &spatialSEIRModel::estimateR0)
    .method("estimateEffectiveR0", &spatialSEIRModel::estimateEffectiveR02)
    .method("estimateEffectiveR0", &spatialSEIRModel::estimateEffectiveR0)
    .method("printAcceptanceRates", &spatialSEIRModel::printAcceptanceRates)
    .method("printOCLSummary", &spatialSEIRModel::printOCLSummary)
    .method("printSamplingParameters", &spatialSEIRModel::printSamplingParameters)
    .method("updateSamplingParameters", &spatialSEIRModel::updateSamplingParameters)
    .method("getGenerationMatrix", &spatialSEIRModel::getGenerationMatrix)
    .method("getIntegratedGenerationMatrix", &spatialSEIRModel::getIntegratedGenerationMatrix)
    .method("standardizeDistanceMatrices", &spatialSEIRModel::standardizeDistanceMatrices)
    .property("S", &spatialSEIRModel::getS, "Susceptible Compartment Matrix")
    .property("E", &spatialSEIRModel::getE, "Exposed Compartment Matrix")
    .property("I", &spatialSEIRModel::getI, "Infectious Compartment Matrix")
    .property("R", &spatialSEIRModel::getR, "Removed Compartment Matrix")
    .property("S0", &spatialSEIRModel::getS0, "Initial Susceptible Compartment Matrix")
    .property("E0", &spatialSEIRModel::getE0, "Initial Exposed Compartment Matrix")
    .property("I0", &spatialSEIRModel::getI0, "Initial Infectious Compartment Matrix")
    .property("R0", &spatialSEIRModel::getR0, "Initial Removed Compartment Matrix")
    .property("S_star", &spatialSEIRModel::getS_star, "Removed to Susceptible Transition Matrix")
    .property("E_star", &spatialSEIRModel::getE_star, "Susceptible to Exposed Transition Matrix")
    .property("I_star", &spatialSEIRModel::getI_star, "Exposed to Infectious Transition Matrix")
    .property("R_star", &spatialSEIRModel::getR_star, "Infectious to Removed Transition Matrix")
    .property("p_se", &spatialSEIRModel::getP_SE, "Exposure Probability Matrix")
    .property("p_ei", &spatialSEIRModel::getP_EI, "E to I Transition Probability")
    .property("p_ir", &spatialSEIRModel::getP_IR, "I to R Transition Probability")
    .property("p_rs", &spatialSEIRModel::getP_RS, "R-S Transition Probability Vector")
    .property("beta", &spatialSEIRModel::getBeta, "Exposure Process Regression Parameters")
    .property("betaP_RS", &spatialSEIRModel::getBetaP_RS, "R-S Transition Process Regression Parameters")
    .property("rho", &spatialSEIRModel::getRho, "Spatial Dependence Term")
    .property("phi", &spatialSEIRModel::getPhi, "Overdispersion Term")
    .property("Y", &spatialSEIRModel::getY, "Raw Data")
    .property("parameterSamplingMode", &spatialSEIRModel::getParameterSamplingMode, 
            &spatialSEIRModel::setParameterSamplingMode, "Type of sampler used for non-compartment parameters.")
    .property("compartmentSamplingMode", &spatialSEIRModel::getCompartmentSamplingMode, 
            &spatialSEIRModel::setCompartmentSamplingMode, "Type of sampler used for disease compartments.")
    .property("debug", &spatialSEIRModel::getDebug, &spatialSEIRModel::setDebug, "Show debug level output?")
    .property("verbose", &spatialSEIRModel::getVerbose, &spatialSEIRModel::setVerbose, "Show verbose level output?")
    .property("useDecorrelation", &spatialSEIRModel::getUseDecorrelation, &spatialSEIRModel::setUseDecorrelation, "Use decorrelation sampling?")
    ;

}

