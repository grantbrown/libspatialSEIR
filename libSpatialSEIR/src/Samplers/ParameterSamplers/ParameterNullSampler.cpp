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
    ParameterNullSampler::ParameterNullSampler() 
    {
    }

    ParameterNullSampler::~ParameterNullSampler()
    {
    }

    int ParameterNullSampler::getSamplerType()
    {
        return(PARAMETER_NULL_SAMPLER);
    }

    void ParameterNullSampler::drawSample()
    {
        // Draw sample code here
    }
}
