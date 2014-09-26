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
    CompartmentBinomialMetropolisSampler::CompartmentBinomialMetropolisSampler(ModelContext* context_,
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

        *context = context_; 
        *compartmentFC = compartmentFC_;
        *compartmentData = compartmentData_;
        *compartmentFrom = compartmentFrom_;
        *compartmentTo = compartmentTo_;
        *probabilityVector = probabilityVector_;
        *probabilityVectorLen = probabilityVectorLen_;
    }

    CompartmentBinomialMetropolisSampler::~CompartmentBinomialMetropolisSampler()
    {
        delete compartmentFC;
        delete compartmentData;
        delete compartmentTo;
        delete compartmentFrom;
        delete probabilityVector;
        delete probabilityVectorLen;
        delete context;
    }

    int CompartmentBinomialMetropolisSampler::getSamplerType()
    {
        return(COMPARTMENT_BINOM_PROPOSAL_METROPOLIS_SAMPLER);
    }

    void CompartmentBinomialMetropolisSampler::drawSample()
    {
        int initAccepted = *((*compartmentFC) -> accepted);
        double initVal;
        int i;
        int x0, x1;
        int nLoc = *((*context) -> S_star -> ncol);
        int nTpt = *((*context) -> S -> nrow);
        int loc;
        int compIdx;
        double proposalNumerator;
        double proposalDenominator;
        double p;
        int n;
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
        for (loc = 0; loc < nLoc; loc++)
        {
            compIdx = loc*nTpt;
            for (i = 0; i < (nTpt); i++)
            { 
                *((*compartmentFC) -> samples) += 1;
                //(*compartmentFC) -> calculateRelevantCompartments(); 
                //(*compartmentFC) -> evalCPU();

                (*compartmentFC) -> calculateRelevantCompartments(loc, i); 
                (*compartmentFC) -> evalCPU(loc, i);
                initVal = (*compartmentFC) -> getValue();
                x0 = (*compartmentData)[compIdx];
                n = (*compartmentFrom)[compIdx];
                p = (*probabilityVector)[compIdx % (*probabilityVectorLen)];
                x1 = (*context) -> random -> binom(n, p);
                proposalNumerator = (*context) -> random -> dbinom(x0, n, p);
                proposalDenominator = (*context) -> random -> dbinom(x1, n, p);
                (*compartmentData)[compIdx] = x1;

                (*compartmentFC) -> calculateRelevantCompartments(loc, i); 
                //(*compartmentFC) -> calculateRelevantCompartments(); 
                //(*compartmentFC) -> evalCPU();
                (*compartmentFC) -> evalCPU(loc, i);
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
                    (*compartmentFC) -> calculateRelevantCompartments(loc, i); 
                    //(*compartmentFC) -> calculateRelevantCompartments(); 

                }

                compIdx++;

            }
        }

        if ((*((*compartmentFC) -> accepted)) == initAccepted)
        {
            // Keep original values
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
