#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_FullConditional.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<LSS_Samplers.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    CompartmentMetropolisSampler::CompartmentMetropolisSampler(ModelContext* context_,
                                                               CompartmentFullConditional* compartmentFC_,
                                                               int* compartmentData_)
    {
        context = new ModelContext*;
        compartmentFC = new CompartmentFullConditional*;
        compartmentData = new int*;

        *context = context_; 
        *compartmentFC = compartmentFC_;
        *compartmentData = compartmentData_;
    }

    CompartmentMetropolisSampler::~CompartmentMetropolisSampler()
    {
        delete compartmentFC;
        delete compartmentData;
        delete context;
    }

    void CompartmentMetropolisSampler::drawSample()
    {
        *((*compartmentFC) -> samples) += 1;
        double initVal;
        double initProposal = 0.0;
        double newProposal = 0.0;
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
            std::cerr << "Compartment sampler starting from value of zero probability.\n";
            throw(-1);
        }
        for (i = 0; i < totalPoints; i++)
        {
            x0 = (*compartmentData)[i];
            x1 = std::floor(((*context) -> random -> normal(x0 + 0.5, sliceWidth)));
            (*compartmentData)[i] = x1;
            newProposal += ((*context) -> random -> dnorm(x1, x0, sliceWidth));
            initProposal += ((*context) -> random -> dnorm(x0, x1, sliceWidth)); 
        }
        (*compartmentFC) -> calculateRelevantCompartments();
        (*compartmentFC) -> evalCPU();

        double newVal = (*compartmentFC) -> getValue();
        double criterion = (newVal - initVal) + (initProposal - newProposal);

        if (std::log((*context) -> random -> uniform()) < criterion)
        {
            // Accept new values
            ((*compartmentFC) -> accepted) += 1;
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
            std::cout << "Impossible value selected.\n";
            throw(-1);
        } 
    }
}
