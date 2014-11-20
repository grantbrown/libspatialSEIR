#include<math.h>
#include<cstring>
#include<vector>
#include<cmath>
#include<algorithm>
#include<LSS_FullConditional.hpp>
#include<LSS_IterationTasks.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>
#include<LSS_FC_Beta.hpp>
#include<LSS_FC_Gamma_EI.hpp>
#include<LSS_Samplers.hpp>
#include<IOProvider.hpp>

namespace SpatialSEIR
{
    PerformHybridSE_EI_UpdateStep::PerformHybridSE_EI_UpdateStep(ModelContext* context_,
                                                                 FC_Gamma_EI* fc_gammaEI_,
                                                                 FC_Beta* fc_beta_,                                                                
                                                                 int iterationCount_)
    {
        context = new ModelContext*;
        iterationCount = new int;
        currentIteration = new int; 

        std::vector<ParameterFullConditional*> fcVector;
        std::vector<double*> params;

        fcVector.push_back(fc_gammaEI_);
        fcVector.push_back(fc_beta_);
        params.push_back(context_ -> gamma_ei);
        params.push_back(context_ -> beta);
        sampler = new ParameterHybridSampler(context_,
                                             fcVector,
                                             params,
                                             HYBRID_SAMPLER_BETA_P_EI);

        *context = context_;
        *iterationCount = iterationCount_;
        *currentIteration = 0;
    }

    int PerformHybridSE_EI_UpdateStep::getTaskType()
    {
        return(LSS_HYBRID_SE_EI_UPDATE_STEP_TASK_TYPE);
    }

    void PerformHybridSE_EI_UpdateStep::executeTask()
    {
    }

    PerformHybridSE_EI_UpdateStep::~PerformHybridSE_EI_UpdateStep()
    {
        delete context;
        delete iterationCount;
        delete currentIteration;
    }
}

