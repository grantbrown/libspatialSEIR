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
    ParameterHybridSampler::ParameterHybridSampler(ModelContext* context_,
                                                   std::vector<ParameterFullConditional*> parameterFullConditionals_,
                                                   std::vector<double*> parameters_,
                                                   int samplerType_) 
    {


        context = new ModelContext*;
        parameterFullConditionals = new std::vector<ParameterFullConditional*>();
        parameters = new std::vector<double*>();
        samplerType = new int;
        *samplerType = samplerType_;
        totalParamSize = new int; *totalParamSize = 0;

        *context = context_;    

        unsigned int i;
        for (i = 0; i < parameterFullConditionals_.size(); i++)
        {
            parameterFullConditionals -> push_back(parameterFullConditionals_[i]);
            parameters -> push_back(parameters_[i]);
            (*totalParamSize) += *(parameterFullConditionals_[i] -> varLen);
        }

        parameterCache = new double[*totalParamSize];
    }

    ParameterHybridSampler::~ParameterHybridSampler()
    {
        delete parameterFullConditionals;
        delete parameters;
        delete samplerType;
        delete[] parameterCache;
        delete totalParamSize;
        delete context;
    }

    int ParameterHybridSampler::getSamplerType()
    {
        return(*samplerType);
    }

    void ParameterHybridSampler::drawSample()
    {
        unsigned int i;
        double sliceWidth;
        double x0, x1;
        int j, k;
        int iters;
        bool success = false;
        long double initVal = 0.0;
        long double newVal = 0.0;
        // Record Current Value
        for (i = 0; i < (parameterFullConditionals -> size()); i++)
        {
            (*parameterFullConditionals)[i] -> evalCPU();
            initVal += (*parameterFullConditionals)[i] -> getValue();
        }

        if (! std::isfinite(initVal))
        {
            lssCout << "Hybrid sampler starting from area of zero probability.\n";
            throw(-1);
        }
        // Back up values
        k = 0;
        for (i = 0; i < (parameterFullConditionals -> size()); i++)
        {
            (*((*parameterFullConditionals)[i] -> samples)) += 1;
            for (j = 0; j < *((*parameterFullConditionals)[i] -> varLen); j++)
            {
                parameterCache[k] = ((*parameters)[i])[j];
                k++;  
            }
        }

        iters = 0;
        while (!success && iters < 1000)
        {
            k = 0;
            newVal = 0.0;
            for (i = 0; i < (parameterFullConditionals -> size()); i++)
            {
                for (j = 0; j < *((*parameterFullConditionals)[i] -> varLen); j++)
                {
                    sliceWidth = *((*parameterFullConditionals)[i] -> sliceWidth);
                    x0 = parameterCache[k];
                    x1 = (((*context) -> random -> normal(x0, sliceWidth)));
                    ((*parameters)[i])[j] = x1;
                    k++;  
                }
            }
            for (i = 0; i < (parameterFullConditionals -> size()); i++)
            {
                (*parameterFullConditionals)[i] -> calculateRelevantCompartments();
                (*parameterFullConditionals)[i] -> evalCPU();
                newVal += (*parameterFullConditionals)[i] -> getValue();
            }
            if (std::log(((*context) -> random -> uniform())) < (newVal - initVal))
            {
                success = true;
                for (i = 0; i < (parameterFullConditionals -> size()); i++)
                {
                    for (j = 0; j < *((*parameterFullConditionals)[i] -> varLen); j++)
                    {
                        ((*parameterFullConditionals)[i] -> accepted)[j]++; 
                    }
                }
            }
            iters ++;
        }
        if (iters >= 1000)
        {
            lssCout << "Hybrid Sampler did Not Update.\n";
            k = 0;
            // Re-set parameter values
            for (i = 0; i < (parameterFullConditionals -> size()); i++)
            {
                (*((*parameterFullConditionals)[i] -> samples)) += 1;
                for (j = 0; j < *((*parameterFullConditionals)[i] -> varLen); j++)
                {
                    ((*parameters)[i])[j] = parameterCache[k];
                    k++;  
                }
                (*parameterFullConditionals)[i] -> calculateRelevantCompartments();
            }
        }
    }
}
