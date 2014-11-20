#include<math.h>
#include<cstring>
#include<vector>
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
    IndexedInitCompartmentMetropolisSampler::IndexedInitCompartmentMetropolisSampler(ModelContext* context_,
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
        *indexLength = (*context) -> indexLength;
        *indexList = (*context) -> indexList;
    }

    IndexedInitCompartmentMetropolisSampler::~IndexedInitCompartmentMetropolisSampler()
    {
        delete initCompartmentFC;
        delete initCompartmentData;
        delete indexLength;
        delete indexList;
        delete context;
    }

    int IndexedInitCompartmentMetropolisSampler::getSamplerType()
    {
        return(INITCOMPARTMENT_IDX_METROPOLIS_SAMPLER);
    }

    void IndexedInitCompartmentMetropolisSampler::drawSample()
    {
        *((*initCompartmentFC) -> samples) += 1;
        double initVal;
        double sliceWidth = *((*initCompartmentFC) -> sliceWidth);
        int i;
        int x0, x1;
        int totalPoints = (*((*context) -> S -> ncol));
        memcpy((*context) -> tmpContainer -> data, *initCompartmentData, totalPoints*sizeof(int));
        (*initCompartmentFC) -> calculateRelevantCompartments(); 
        (*initCompartmentFC) -> evalCPU();
        initVal = (*initCompartmentFC) -> getValue();
        if (! std::isfinite(initVal))
        {
            lssCout << "Init compartment sampler starting from value of zero probability.\n";
            throw(-1);
        }
        for (i = 0; i < **indexLength; i++)
        {
            x0 = (*initCompartmentData)[(*indexList)[i]];
            x1 = std::floor(((*context) -> random -> normal(x0 + 0.5, sliceWidth)));
            (*initCompartmentData)[(*indexList)[i]] = x1;
        }
        (*initCompartmentFC) -> calculateRelevantCompartments(); 
        (*initCompartmentFC) -> evalCPU();
        double newVal = (*initCompartmentFC) -> getValue();
        double criterion = (newVal - initVal);

        if (std::log((*context) -> random -> uniform()) < criterion)
        {
            // Accept new values
            *((*initCompartmentFC) -> accepted) += 1;
        }
        else
        {
            // Keep original values
            memcpy(*initCompartmentData, (*context) -> tmpContainer -> data, totalPoints*sizeof(int));
            (*initCompartmentFC) -> calculateRelevantCompartments(); 
            (*initCompartmentFC) -> setValue(initVal); 
        }
        if (! std::isfinite((*initCompartmentFC) -> getValue()))
        {
            lssCout << "Impossible value selected.\n";
            throw(-1);
        } 


    }
}
