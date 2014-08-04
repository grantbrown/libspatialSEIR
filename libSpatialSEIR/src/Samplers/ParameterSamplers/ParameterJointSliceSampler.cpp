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
    ParameterJointSliceSampler::ParameterJointSliceSampler(ModelContext* context_,
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

    ParameterJointSliceSampler::~ParameterJointSliceSampler()
    {
        delete paramFC;
        delete param;
        delete context;
    }

    int ParameterJointSliceSampler::getSamplerType()
    {
        return(PARAMETER_JOINT_SLICE_SAMPLER);
    }

    void ParameterJointSliceSampler::drawSample()
    {
        lssCout << "Joint slice sampling not yet implemented\n";
        throw(-1); 
    }
}
