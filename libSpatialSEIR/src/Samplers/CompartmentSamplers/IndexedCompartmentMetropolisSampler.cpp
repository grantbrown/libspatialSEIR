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
#include<IOProvider.hpp>
namespace SpatialSEIR
{
    IndexedCompartmentMetropolisSampler::IndexedCompartmentMetropolisSampler(ModelContext* context_,
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

    IndexedCompartmentMetropolisSampler::~IndexedCompartmentMetropolisSampler()
    {
        delete compartmentFC;
        delete compartmentData;
        delete indexLength;
        delete indexList;
        delete context;
    }

    int IndexedCompartmentMetropolisSampler::getSamplerType()
    {
        return(COMPARTMENT_IDX_METROPOLIS_SAMPLER);
    }

    void IndexedCompartmentMetropolisSampler::drawSample()
    {
        *((*compartmentFC) -> samples) += 1;
        double initVal;
        double sliceWidth = *((*compartmentFC) -> sliceWidth);
        int i;
        int x0, x1;
        int totalPoints = *((*context) -> S -> nrow)*(*((*context) -> S -> ncol));
        memcpy((*context) -> tmpContainer -> data, *compartmentData, totalPoints*sizeof(int));
        (*compartmentFC) -> calculateRelevantCompartments(); 
        (*compartmentFC) -> evalCPU();
        initVal = (*compartmentFC) -> getValue();
        if (! std::isfinite(initVal))
        {
            lssCout << "Compartment sampler starting from value of zero probability.\n";
            throw(-1);
        }
        for (i = 0; i < (**indexLength); i++)
        {
            x0 = (*compartmentData)[(*indexList)[i]];
            x1 = std::floor(((*context) -> random -> normal(x0 + 0.5, sliceWidth)));
            (*compartmentData)[(*indexList)[i]] = x1;
        }
        (*compartmentFC) -> calculateRelevantCompartments();
        (*compartmentFC) -> evalCPU();

        double newVal = (*compartmentFC) -> getValue();
        double criterion = (newVal - initVal);

        if (std::log((*context) -> random -> uniform()) < criterion)
        {
            // Accept new values
            *((*compartmentFC) -> accepted) += 1;
        }
        else
        {
            // Keep original values
            memcpy(*compartmentData, (*context) -> tmpContainer -> data, totalPoints*sizeof(int));
            (*compartmentFC) -> calculateRelevantCompartments();
            (*compartmentFC) -> setValue(initVal);
        }
        if (! std::isfinite((*compartmentFC) -> getValue()))
        {
            lssCout << "Impossible value selected.\n";
            throw(-1);
        } 
    }
}
