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
#include<IOProvider.hpp>
namespace SpatialSEIR
{
    ParameterHybridSampler::ParameterHybridSampler(ModelContext* context_,
                                                   std::vector<ParameterFullConditional*> parameterFullConditionals_,
                                                   std::vector<double*> parameters_,
                                                   int samplerType_) 
    {

        context = new ModelContext*;
        parameterFullConditionals = new std::vector<ParameterFullConditional*>;
        parameters = new std::vector<double*>;
        samplerType = new int;
        *samplerType = samplerType_;

        *context = context_;    
        unsigned int i;
        for (i = 0; i < parameterFullConditionals_.size(); i++)
        {
            parameterFullConditionals -> push_back(parameterFullConditionals_[i]);
            parameters -> push_back(parameters_[i]);
        }
    }

    ParameterHybridSampler::~ParameterHybridSampler()
    {
        delete parameterFullConditionals;
        delete parameters;
        delete samplerType;
        delete context;
    }

    int ParameterHybridSampler::getSamplerType()
    {
        return(*samplerType);
    }

    void ParameterHybridSampler::drawSample()
    {
        // Draw sample code here
    }
}
