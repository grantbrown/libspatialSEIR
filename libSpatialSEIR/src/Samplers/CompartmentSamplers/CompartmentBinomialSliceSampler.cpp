#include<math.h>
#include<cstring>
#include<vector>
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
    CompartmentBinomialSliceSampler::CompartmentBinomialSliceSampler(ModelContext* context_,
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

    CompartmentBinomialSliceSampler::~CompartmentBinomialSliceSampler()
    {
        delete compartmentFC;
        delete compartmentData;
        delete compartmentTo;
        delete compartmentFrom;
        delete probabilityVector;
        delete probabilityVectorLen;
        delete context;
    }

    int CompartmentBinomialSliceSampler::getSamplerType()
    {
        return(COMPARTMENT_BINOM_PROPOSAL_SLICE_SAMPLER);
    }

    void CompartmentBinomialSliceSampler::drawSample()
    {
        double initVal;
        int i;
        //int x0, x1;
        int nLoc = *((*context) -> S_star -> ncol);
        int nTpt = *((*context) -> S -> nrow);
        int loc;
        int compIdx;
        double p;
        int n;
        double l, r, y, x0, x1;

        (*compartmentFC) -> calculateRelevantCompartments(); 
        (*compartmentFC) -> evalCPU();
        initVal = (*compartmentFC) -> getValue();
        if (! std::isfinite(initVal))
        {
            lssCout << "Compartment sampler starting from value of zero probability.\n";
            throw(-1);
        }

        for (loc = 0; loc < nLoc; loc++)
        {
            compIdx = loc*nTpt;
            for (i = 0; i < (nTpt); i++)
            { 
                (*compartmentFC) -> calculateRelevantCompartments(loc, i); 
                (*compartmentFC) -> evalCPU(loc, i);
                initVal = (*compartmentFC) -> getValue();
                y = initVal - ((*context) -> random -> gamma());

                x0 = (*compartmentData)[compIdx];
                n = (*compartmentFrom)[compIdx];
                p = (*probabilityVector)[compIdx % (*probabilityVectorLen)];
                x1 = (*context) -> random -> binom(n, p);
                (*compartmentData)[compIdx] = x1;

                (*compartmentFC) -> calculateRelevantCompartments(loc, i); 
                (*compartmentFC) -> evalCPU(loc, i);

                while (y >= (*compartmentFC) -> getValue())
                {
                    l = std::min(x0, x1);
                    r = std::max(x0, x1) + 1;
                    x1 = (((*context) -> random -> uniform())*(r-l) + l);
                    (*compartmentData)[compIdx] = std::floor(x1);
                    (*compartmentFC) -> calculateRelevantCompartments(loc, i); 
                    (*compartmentFC) -> evalCPU(loc, i);
                    if (x1 >= x0){r=x1;}
                    else{l=x1;}
                }
                compIdx++;
            }
        }

        if (! std::isfinite((*compartmentFC) -> getValue()))
        {
            lssCout << "Impossible value selected.\n";
            throw(-1);
        } 
    }
}
