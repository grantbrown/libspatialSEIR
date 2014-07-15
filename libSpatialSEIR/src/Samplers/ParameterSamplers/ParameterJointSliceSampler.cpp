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
        std::cout << "Joint slice sampling not yet implemented\n";
        throw(-1); 
    }
}