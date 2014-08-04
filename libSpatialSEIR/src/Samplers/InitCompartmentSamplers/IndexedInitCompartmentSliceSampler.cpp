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
    IndexedInitCompartmentSliceSampler::IndexedInitCompartmentSliceSampler(ModelContext* context_,
                                                               InitCompartmentFullConditional* initCompartmentFC_,
                                                               int* initCompartmentData_)
    {
        context = new ModelContext*;
        initCompartmentFC = new InitCompartmentFullConditional*;
        initCompartmentData = new int*;
        indexList = new int*;
        indexLength = new int*;

        *context = context_; 
        *initCompartmentFC = initCompartmentFC_;
        *initCompartmentData = initCompartmentData_;
        *indexLength = (*context) -> indexList;
        *indexList = (*context) -> indexLength;
    }

    IndexedInitCompartmentSliceSampler::~IndexedInitCompartmentSliceSampler()
    {
        delete initCompartmentFC;
        delete initCompartmentData;
        delete indexLength;
        delete indexList;
        delete context;
    }

    int IndexedInitCompartmentSliceSampler::getSamplerType()
    {
        return(INITCOMPARTMENT_IDX_SLICE_SAMPLER);
    }

    void IndexedInitCompartmentSliceSampler::drawSample()
    {
        lssCout << "Blocked slice sampling not yet implemented.\n"; 
        throw(-1);
    }
}
