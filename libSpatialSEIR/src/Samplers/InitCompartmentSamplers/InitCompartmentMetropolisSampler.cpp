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

    int InitCompartmentMetropolisSampler::getSamplerType()
    {
        return(INITCOMPARTMENT_METROPOLIS_SAMPLER);
    }

    void InitCompartmentMetropolisSampler::drawSample()
    {



        *((*initCompartmentFC) -> samples) += 1;
        double initVal;
        double initProposal = 0.0;
        double newProposal = 0.0;
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
            std::cerr << "Compartment sampler starting from value of zero probability.\n";
            throw(-1);
        }
        for (i = 0; i < totalPoints; i++)
        {
            x0 = (*initCompartmentData)[i];
            x1 = std::floor(((*context) -> random -> normal(x0 + 0.5, sliceWidth)));
            (*initCompartmentData)[i] = x1;
            newProposal += ((*context) -> random -> dnorm(x1, x0 + 0.5, sliceWidth));
            initProposal += ((*context) -> random -> dnorm(x0, x1 + 0.5, sliceWidth)); 
        }
        (*initCompartmentFC) -> calculateRelevantCompartments(); 
        (*initCompartmentFC) -> evalCPU();
        double newVal = (*initCompartmentFC) -> getValue();
        double criterion = (newVal - initVal) + (initProposal - newProposal);

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
            std::cout << "Impossible value selected.\n";
            throw(-1);
        } 

        
    }
}
