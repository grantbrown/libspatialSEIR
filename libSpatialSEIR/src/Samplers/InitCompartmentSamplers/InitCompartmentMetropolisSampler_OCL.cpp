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
    InitCompartmentMetropolisSampler_OCL::InitCompartmentMetropolisSampler_OCL(ModelContext* context_,
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

    InitCompartmentMetropolisSampler_OCL::~InitCompartmentMetropolisSampler_OCL()
    {
        delete initCompartmentFC;
        delete initCompartmentData;
        delete context;
    }

    int InitCompartmentMetropolisSampler_OCL::getSamplerType()
    {
        return(INITCOMPARTMENT_METROPOLIS_SAMPLER_OCL);
    }

    void InitCompartmentMetropolisSampler_OCL::drawSample()
    {



        *((*initCompartmentFC) -> samples) += 1;
        double initVal;
        double sliceWidth = *((*initCompartmentFC) -> sliceWidth);
        int i;
        int x0, x1;
        int totalPoints = (*((*context) -> S -> ncol));
        memcpy((*context) -> tmpContainer -> data, *initCompartmentData, totalPoints*sizeof(int));
        (*initCompartmentFC) -> calculateRelevantCompartments_OCL(); 
        (*initCompartmentFC) -> evalOCL();
        initVal = (*initCompartmentFC) -> getValue();
        if (! std::isfinite(initVal))
        {
            lssCout << "Init compartment sampler starting from value of zero probability.\n";
            throw(-1);
        }
        for (i = 0; i < totalPoints; i++)
        {
            x0 = (*initCompartmentData)[i];
            x1 = std::floor(((*context) -> random -> normal(x0 + 0.5, sliceWidth)));
            (*initCompartmentData)[i] = x1;
        }
        (*initCompartmentFC) -> calculateRelevantCompartments_OCL(); 
        (*initCompartmentFC) -> evalOCL();
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
            (*initCompartmentFC) -> calculateRelevantCompartments_OCL(); 
            (*initCompartmentFC) -> setValue(initVal); 
        }
        if (! std::isfinite((*initCompartmentFC) -> getValue()))
        {
            lssCout << "Impossible value selected.\n";
            throw(-1);
        } 

        
    }
}
