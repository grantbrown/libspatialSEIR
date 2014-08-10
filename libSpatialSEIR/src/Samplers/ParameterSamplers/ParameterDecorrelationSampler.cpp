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
#include<cblas.h>

namespace SpatialSEIR
{
    ParameterDecorrelationSampler::ParameterDecorrelationSampler(ModelContext* context_,
                                                                       ParameterFullConditional* paramFC_,
                                                                       double* param_,
                                                                       CovariateMatrix* proposalMatrix_) 
    {
        context = new ModelContext*;
        paramFC = new ParameterFullConditional*;
        param = new double*;
        proposalCache = new double[*((*paramFC) -> varLen)];
        proposalMatrix = new CovariateMatrix*;

        *context = context_;    
        *paramFC = paramFC_;
        *param = param_;
        *proposalMatrix = proposalMatrix_;
    }

    ParameterDecorrelationSampler::~ParameterDecorrelationSampler()
    {
        delete paramFC;
        delete param;
        delete context;
        delete proposalMatrix;
        delete[] proposalCache;
    }

    int ParameterDecorrelationSampler::getSamplerType()
    {
        return(PARAMETER_DECORR_SAMPLER);
    }

    void ParameterDecorrelationSampler::drawSample()
    {
        *((*paramFC) -> samples) += 1;
        double initVal;
        double sliceWidth = *((*paramFC) -> sliceWidth);
        int i;
        int totalPoints = *((*paramFC) -> varLen);
        memcpy((*context) -> compartmentCache, *param, totalPoints*sizeof(double));
        (*paramFC) -> calculateRelevantCompartments(); 
        (*paramFC) -> evalCPU();
        initVal = (*paramFC) -> getValue();

        if (! std::isfinite(initVal))
        {
            lssCout << "Parameter sampler starting from value of zero probability.\n";
            throw(-1);
        }

        for (i = 0; i < totalPoints; i++)
        {
            proposalCache[i] = (((*context) -> random -> normal(0, 1))); 
        }

        cblas_dgemv(CblasColMajor, 
                    CblasNoTrans,
                    totalPoints,
                    totalPoints,
                    1.0,
                    ((*proposalMatrix) -> decorrelationProjectionMatrix),
                    totalPoints,
                    proposalCache,
                    1,
                    0.0,
                    proposalCache,
                    1);

        // End updates.
        double proposalCacheSize = 0.;
        for (i = 0; i < totalPoints; i++)
        {
            proposalCacheSize += (proposalCache[i])*(proposalCache[i]);
        }
        proposalCacheSize = std::sqrt(proposalCacheSize);
        for (i = 0; i < totalPoints; i++)
        {
            proposalCache[i] /= (proposalCacheSize/sliceWidth);
            (*param)[i] += proposalCache[i];
        }

        (*paramFC) -> calculateRelevantCompartments(); 
        (*paramFC) -> evalCPU();
        double newVal = (*paramFC) -> getValue();
        double criterion = (newVal - initVal);

        if (std::log((*context) -> random -> uniform()) < criterion)
        {
            // Accept new values
            for (i = 0; i < totalPoints; i++)
            {
                ((*paramFC) -> accepted)[i] += 1;
            }
        }
        else
        {
            // Keep original values
            memcpy(*param, (*context) -> compartmentCache, totalPoints*sizeof(double));
            (*paramFC) -> calculateRelevantCompartments(); 
            (*paramFC) -> setValue(initVal); 
        }

        if (! std::isfinite((*paramFC) -> getValue()))
        {
            lssCout << "Impossible value selected.\n";
            throw(-1);
        }
    }
}
