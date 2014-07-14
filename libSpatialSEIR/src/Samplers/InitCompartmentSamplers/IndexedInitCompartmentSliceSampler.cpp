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

    IndexedInitCompartmentSliceSampler::IndexedInitCompartmentSliceSampler(ModelContext* context_,
                                                               InitCompartmentFullConditional* initCompartmentFC_,
                                                               int* initCompartmentData_,
                                                               int* indexList_,
                                                               int indexLength_)
    {
        context = new ModelContext*;
        initCompartmentFC = new InitCompartmentFullConditional*;
        initCompartmentData = new int*;
        indexList = new int*;
        indexLength = new int;

        *context = context_; 
        *initCompartmentFC = initCompartmentFC_;
        *initCompartmentData = initCompartmentData_;
        *indexLength = indexLength_;
        *indexList = indexList_;
    }

    IndexedInitCompartmentSliceSampler::~IndexedInitCompartmentSliceSampler()
    {
        delete initCompartmentFC;
        delete initCompartmentData;
        delete indexLength;
        delete indexList;
        delete context;
    }

    void IndexedInitCompartmentSliceSampler::drawSample()
    {
        
    }
}
