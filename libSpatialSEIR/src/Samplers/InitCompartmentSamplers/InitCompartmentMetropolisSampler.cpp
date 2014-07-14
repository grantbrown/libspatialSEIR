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

    InitCompartmentMetropolisSampler::InitCompartmentMetropolisSampler(ModelContext* context_,
                                                               InitCompartmentFullConditional* initCompartmentFC_,
                                                               int* initCompartmentData_)
    {
        context = new ModelContext*;
        initCompartmentFC = new InitCompartmentFullConditional*;
        initCompartmentData = new int*;

        *context = context_; 
        *initCompartmentFC = initCompartmentFC_;
        *initCompartmentData = initCompartmentData_;
    }

    InitCompartmentMetropolisSampler::~InitCompartmentMetropolisSampler()
    {
        delete initCompartmentFC;
        delete initCompartmentData;
        delete context;
    }

    void InitCompartmentMetropolisSampler::drawSample()
    {
        
    }
}
