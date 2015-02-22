#include<math.h>
#include<cstring>
#include<vector>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FullConditional.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>
#include<IOProvider.hpp>
namespace SpatialSEIR
{
    ParameterSingleMetropolisSampler::ParameterSingleMetropolisSampler(ModelContext* context_,
                                                                       ParameterFullConditional* paramFC_,
                                                                       double* param_) 
    {
        context = new ModelContext*;
        paramFC = new ParameterFullConditional*;
        param = new double*;

        *context = context_;    
        *paramFC = paramFC_;
        *param = param_;
    }

    ParameterSingleMetropolisSampler::~ParameterSingleMetropolisSampler()
    {
        delete paramFC;
        delete param;
        delete context;
    }

    int ParameterSingleMetropolisSampler::getSamplerType()
    {
        return(PARAMETER_SINGLE_METROPOLIS_SAMPLER);
    }

    void ParameterSingleMetropolisSampler::drawSample()
    {
        // Declare required variables
        int i;
        int varLen = *((*paramFC) -> varLen);
        double* width = ((*paramFC) -> sliceWidth);
        double x0,x1;
        double initVal, newVal;

        (*((*paramFC) -> samples)) += 1;
        // Update the relevant CompartmentalModelMatrix instances
        (*paramFC) -> calculateRelevantCompartments();

        // Set the "value" attribute appropriately
        (*paramFC) -> evalCPU();
        initVal = ((*paramFC) -> getValue());
        if (! std::isfinite(initVal))
        {
            lssCout << "Parameter sampler starting from value of zero probability.\n";
            throw(-1);
        }
  
        // Main loop: 
        for (i = 0; i < varLen; i++)
        { 
            x0 = (*param)[i];
            (*paramFC) -> calculateRelevantCompartments(); 
            (*paramFC) -> evalCPU();
            initVal = ((*paramFC)->getValue());

            x1 = ((*context) -> random -> normal(x0, width[i]));
            (*param)[i] = x1;
            (*paramFC) -> calculateRelevantCompartments();
            (*paramFC) -> evalCPU();
            newVal = ((*paramFC)->getValue());

            if (std::log(((*context) -> random -> uniform())) < ((newVal - initVal)))
            {
                // Accept the new value. 
                ((*paramFC) -> accepted)[i]+=1;
            }
            else
            {
                // Keep original value
                (*param)[i] = x0;
                (*paramFC) -> calculateRelevantCompartments();
                (*paramFC) -> setValue(initVal);
            }
        }
    }
}
