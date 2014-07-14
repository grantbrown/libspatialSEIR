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

    ParameterDecorrelationSampler::ParameterDecorrelationSampler(ModelContext* context_,
                                                                       ParameterFullConditional* paramFC_,
                                                                       double* param_,
                                                                       CovariateMatrix* proposalMatrix_) 
    {
        context = new ModelContext*;
        paramFC = new ParameterFullConditional*;
        param = new double*;
        proposalMatrix = new CovariateMatrix*;

        *context = context_;    
        *paramFC = paramFC_;
        *param = param_;
        *proposalMatrix = proposalMatrix_;
    }

    ParameterDecorrelationSampler::~ParameterDecorrelationSampler()
    {
        delete paramFC;
        delete param;
        delete context;
        delete proposalMatrix;
    }

    void ParameterDecorrelationSampler::drawSample()
    {
        
    }
}
