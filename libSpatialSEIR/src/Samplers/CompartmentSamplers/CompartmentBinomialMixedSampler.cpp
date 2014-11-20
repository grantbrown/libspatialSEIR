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
#include<IOProvider.hpp>

namespace SpatialSEIR
{
    CompartmentBinomialMixedSampler::CompartmentBinomialMixedSampler(ModelContext* context_,
                                                               CompartmentFullConditional* compartmentFC_,
                                                               int* compartmentData_,
                                                               int* compartmentFrom_,
                                                               int* compartmentTo_,
                                                               double* probabilityVector_,
                                                               int probabilityVectorLen_)
    {
        context = new ModelContext*;
        compartmentFC = new CompartmentFullConditional*;
        compartmentData = new int*;
        compartmentTo = new int*;
        compartmentFrom = new int*;
        probabilityVector = new double*;
        probabilityVectorLen = new int;
        indexList = new int*;
        indexLength = new int*;
        currentSamplerType = new int;

        *context = context_; 
        *compartmentFC = compartmentFC_;
        *compartmentData = compartmentData_;
        *compartmentFrom = compartmentFrom_;
        *compartmentTo = compartmentTo_;
        *probabilityVector = probabilityVector_;
        *probabilityVectorLen = probabilityVectorLen_;
        *indexLength = (*context) -> indexLength;
        *indexList = (*context) -> indexList;
        *currentSamplerType = 0;
    }

    CompartmentBinomialMixedSampler::~CompartmentBinomialMixedSampler()
    {
        delete compartmentFC;
        delete compartmentData;
        delete compartmentTo;
        delete compartmentFrom;
        delete probabilityVector;
        delete probabilityVectorLen;
        delete context;
        delete indexLength;
        delete indexList;
        delete currentSamplerType;
    }

    int CompartmentBinomialMixedSampler::getSamplerType()
    {
        return(COMPARTMENT_BINOM_MIXED_SAMPLER);
    }

    void CompartmentBinomialMixedSampler::drawSample()
    {
        if (*currentSamplerType == 0)
        {
            *currentSamplerType = 1;
            this -> drawNormalProposalSample();
            return;
        }
        *currentSamplerType = 0;
        this -> drawBinomialProposalSample(); 
    }

    void CompartmentBinomialMixedSampler::drawNormalProposalSample()
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
        for (i = 0; i < totalPoints; i++)
        {
            x0 = (*compartmentData)[i];
            x1 = std::floor(((*context) -> random -> normal(x0 + 0.5, sliceWidth)));
            (*compartmentData)[i] = x1;
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

    void CompartmentBinomialMixedSampler::drawBinomialProposalSample()
    {
        int initAccepted = *((*compartmentFC) -> accepted);
        int nTpt = *((*context) -> S -> nrow);
        int loc, tpt, idxLen, compIdx, x0, x1, i, n;
        double initVal, proposalNumerator, proposalDenominator, p;
        //memcpy((*context) -> tmpContainer -> data, &((*compartmentData)[compIdx]), nTpt*sizeof(int));
        //memcpy(((*context) -> tmpContainer -> data) + nTpt*sizeof(int), &((*compartmentFrom)[compIdx]), nTpt*sizeof(int));

        (*compartmentFC) -> calculateRelevantCompartments(); 
        (*compartmentFC) -> evalCPU();
        initVal = (*compartmentFC) -> getValue();
        if (! std::isfinite(initVal))
        {
            lssCout << "Compartment sampler starting from value of zero probability.\n";
            throw(-1);
        }

        proposalNumerator = 0.0;
        proposalDenominator = 0.0;
        idxLen = **indexLength;
        for (i = 0; i < idxLen; i++)
        {
            compIdx = (*indexList)[i];
            loc = compIdx/(nTpt); 
            tpt = compIdx - loc*nTpt;

            *((*compartmentFC) -> samples) += 1;
            (*compartmentFC) -> calculateRelevantCompartments(loc, tpt); 
            (*compartmentFC) -> evalCPU(loc, tpt);
            initVal = (*compartmentFC) -> getValue();

            x0 = (*compartmentData)[compIdx];
            n = (*compartmentFrom)[compIdx];
            p = (*probabilityVector)[compIdx % (*probabilityVectorLen)];
            x1 = (*context) -> random -> binom(n, p);
            proposalNumerator = (*context) -> random -> dbinom(x0, n, p);
            proposalDenominator = (*context) -> random -> dbinom(x1, n, p);
            (*compartmentData)[compIdx] = x1;

            (*compartmentFC) -> calculateRelevantCompartments(loc, tpt); 
            (*compartmentFC) -> evalCPU(loc, tpt);
            double newVal = (*compartmentFC) -> getValue();
            double criterion = (newVal - initVal) + (proposalNumerator - proposalDenominator);


            if (std::log((*context) -> random -> uniform()) < criterion)
            {
                // Accept new values
                *((*compartmentFC) -> accepted) += 1;
            }
            else 
            {
                (*compartmentData)[compIdx] = x0;
                (*compartmentFC) -> calculateRelevantCompartments(loc, tpt); 
                (*compartmentFC) -> setValue(initVal);
            }
        }

        if ((*((*compartmentFC) -> accepted)) == initAccepted)
        {
            // Keep original values
            (*compartmentFC) -> calculateRelevantCompartments(); 
            (*compartmentFC) -> evalCPU();             
        }

        if (! std::isfinite((*compartmentFC) -> getValue()))
        {
            lssCout << "Impossible value selected: " << (*compartmentFC) -> getValue() << "\n";
            throw(-1);
        } 
    }
}
