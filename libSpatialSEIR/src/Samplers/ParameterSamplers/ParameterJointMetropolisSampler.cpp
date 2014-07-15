#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FullConditional.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    ParameterJointMetropolisSampler::ParameterJointMetropolisSampler(ModelContext* context_,
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

    ParameterJointMetropolisSampler::~ParameterJointMetropolisSampler()
    {
        delete paramFC;
        delete param;
        delete context;
    }

    int ParameterJointMetropolisSampler::getSamplerType()
    {
        return(PARAMETER_JOINT_METROPOLIS_SAMPLER);
    }

    void ParameterJointMetropolisSampler::drawSample()
    {
        *((*paramFC) -> samples) += 1;
        double initVal;
        double sliceWidth = *((*paramFC) -> sliceWidth);
        int i;
        double x0, x1;
        int totalPoints = *((*paramFC) -> varLen);
        memcpy((*context) -> compartmentCache, *param, totalPoints*sizeof(double));
        (*paramFC) -> calculateRelevantCompartments(); 
        (*paramFC) -> evalCPU();
        initVal = (*paramFC) -> getValue();
        if (! std::isfinite(initVal))
        {
            std::cerr << "Compartment sampler starting from value of zero probability.\n";
            throw(-1);
        }
        for (i = 0; i < totalPoints; i++)
        {
            x0 = (*param)[i];
            x1 = std::floor(((*context) -> random -> normal(x0, sliceWidth)));
            (*param)[i] = x1;
        }
        (*paramFC) -> calculateRelevantCompartments(); 
        (*paramFC) -> evalCPU();
        double newVal = (*paramFC) -> getValue();
        double criterion = (newVal - initVal);

        if (std::log((*context) -> random -> uniform()) < criterion)
        {
            // Accept new values
            ((*paramFC) -> accepted) += 1;
        }
        else
        {
            // Keep original values
            memcpy(*param, (*context) -> compartmentCache, totalPoints*sizeof(double));
            (*paramFC) -> calculateRelevantCompartments(); 
            (*paramFC) -> setValue(initVal); 
        }
        if (! std::isfinite((*paramFC) -> getValue()))
        {
            std::cout << "Impossible value selected.\n";
            throw(-1);
        } 
        
    }
}
