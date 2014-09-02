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
                                                               int probabilityVectorLength_)
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
        *probabilityVEctorLen = probabilityVectorLen_;
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
        *((*compartmentFC) -> samples) += 1;
        int initAccepted = *((*compartmentFC) -> accepted);
        double initVal;
        int i;
        int itrs = 0;
        int x0, x1;
        int nLoc = *((*context) -> S_star -> ncol);
        int loc = std::floor((*context) -> random -> uniform(0,nLoc))
        int totalPoints = *((*context) -> S -> nrow);
        int compIdx;
        double proposalNumerator;
        double proposalDenominator;
        double p;
        int n;
        memcpy((*context) -> tmpContainer -> data, *compartmentData, totalPoints*sizeof(int));
        (*compartmentFC) -> calculateRelevantCompartments(); 
        (*compartmentFC) -> evalCPU();
        initVal = (*compartmentFC) -> getValue();

        if (! std::isfinite(initVal))
        {
            lssCout << "Compartment sampler starting from value of zero probability.\n";
            throw(-1);
        }
        do{ 
            proposalNumerator = 0.0;
            proposalDenominator = 0.0;
            for (i = 0; i < totalPoints; i++)
            {
                compIdx = loc*totalPoints + i;
                x0 = (*compartmentData)[compIdx];
                n = (*compartmentFrom)[compIdx];
                p = (*probabilityVector)[compIdx % * probabilityVectorLen];
                //x1 = std::floor(((*context) -> random -> normal(x0 + 0.5, sliceWidth)));
                x1 = (*context) -> random -> binom(n, p);
                proposalNumerator += (*context) -> random -> dbinom(x0, n, p);
                proposalDenominator += (*context) -> random -> dbinom(x1, n, p);
                (*compartmentData)[compIdx] = x1;
            }
            (*compartmentFC) -> calculateRelevantCompartments(); 
            (*compartmentFC) -> evalCPU();
            double newVal = (*compartmentFC) -> getValue();
            double criterion = (newVal - initVal) + (proposalNumerator - proposalDenominator);

            if (std::log((*context) -> random -> uniform()) < criterion)
            {
                // Accept new values
                *((*compartmentFC) -> accepted) += 1;
            }
            itrs++;
        } while(itrs < 1000 && (*((*compartmentFC) -> accepted)) == initAccepted);


        if ((*((*compartmentFC) -> accepted)) == initAccepted)
        {
            // Keep original values
            lssCout << "Compartment sampler did not update.\n";
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
