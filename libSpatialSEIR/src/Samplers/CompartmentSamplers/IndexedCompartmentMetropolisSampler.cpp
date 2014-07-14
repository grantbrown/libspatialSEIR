#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_FullConditional.hpp>
#include<LSS_Samplers.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    IndexedCompartmentMetropolisSampler::IndexedCompartmentMetropolisSampler(ModelContext* context_,
                                                               CompartmentFullConditional* compartmentFC_,
                                                               int* compartmentData_,
                                                               int* indexList_,
                                                               int indexLength_)
    {
        context = new ModelContext*;
        compartmentFC = new CompartmentFullConditional*;
        compartmentData = new int*;
        indexList = new int*;
        indexLength = new int;

        *context = context_; 
        *compartmentFC = compartmentFC_;
        *compartmentData = compartmentData_;
        *indexLength = indexLength_;
        *indexList = indexList_;
    }

    IndexedCompartmentMetropolisSampler::~IndexedCompartmentMetropolisSampler()
    {
        delete compartmentFC;
        delete compartmentData;
        delete indexLength;
        delete indexList;
        delete context;
    }

    void IndexedCompartmentMetropolisSampler::drawSample()
    {
        
    }
}
