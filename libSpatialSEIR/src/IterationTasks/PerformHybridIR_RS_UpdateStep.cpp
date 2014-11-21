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
#include<LSS_FC_Beta_P_RS.hpp>
#include<LSS_FC_Gamma_IR.hpp>
#include<LSS_Samplers.hpp>
#include<IOProvider.hpp>

namespace SpatialSEIR
{
    PerformHybridIR_RS_UpdateStep::PerformHybridIR_RS_UpdateStep(ModelContext* context_,
                                                                 FC_Gamma_IR* fc_gammaIR_,
                                                                 FC_Beta_P_RS* fc_betaPrs_, 
                                                                 int iterationCount_)
    {
        context = new ModelContext*;
        iterationCount = new int;
        currentIteration = new int; 

        std::vector<ParameterFullConditional*> fcVector;
        std::vector<double*> params;

        fcVector.push_back(fc_gammaIR_);
        fcVector.push_back(fc_betaPrs_);
        params.push_back(context_ -> gamma_ir);
        params.push_back(context_ -> betaPrs);
        sampler = new ParameterHybridSampler(context_,
                                             fcVector,
                                             params,
                                             HYBRID_SAMPLER_BETA_RS_P_IR);

        *context = context_;
        *iterationCount = iterationCount_;
        *currentIteration = 0;
    }

    int PerformHybridIR_RS_UpdateStep::getTaskType()
    {
        return(LSS_HYBRID_IR_RS_UPDATE_STEP_TASK_TYPE);
    }

    void PerformHybridIR_RS_UpdateStep::executeTask()
    {
        int hybridEIStride = ((*context) -> config -> performHybridStep);
        if ((*currentIteration) % hybridEIStride != 0)
        {
            *currentIteration += 1;
            return;
        }
        *currentIteration = 1;
        
        (sampler -> drawSample());
    }

    PerformHybridIR_RS_UpdateStep::~PerformHybridIR_RS_UpdateStep()
    {
        delete context;
        delete iterationCount;
        delete currentIteration;
    }
}

