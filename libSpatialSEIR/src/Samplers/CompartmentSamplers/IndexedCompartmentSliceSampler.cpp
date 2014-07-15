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

    IndexedCompartmentSliceSampler::IndexedCompartmentSliceSampler(ModelContext* context_,
                                                               CompartmentFullConditional* compartmentFC_,
                                                               int* compartmentData_)
    {
        context = new ModelContext*;
        compartmentFC = new CompartmentFullConditional*;
        compartmentData = new int*;
        indexList = new int*;
        indexLength = new int*;

        *context = context_; 
        *compartmentFC = compartmentFC_;
        *compartmentData = compartmentData_;
        *indexLength = (*context) -> indexLength;
        *indexList = (*context) -> indexList;
    }

    IndexedCompartmentSliceSampler::~IndexedCompartmentSliceSampler()
    {
        delete compartmentFC;
        delete compartmentData;
        delete indexLength;
        delete indexList;
        delete context;
    }

    int IndexedCompartmentSliceSampler::getSamplerType()
    {
        return(COMPARTMENT_IDX_SLICE_SAMPLER);
    }

    void IndexedCompartmentSliceSampler::drawSample()
    {
        std::cerr << "Block slice sampling not yet implemented.\n";
        throw(-1); 
    }
}
