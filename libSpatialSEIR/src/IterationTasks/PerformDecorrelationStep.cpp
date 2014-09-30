#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
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
#include<LSS_FC_Beta_P_RS.hpp>
#include<LSS_Samplers.hpp>
#include<IOProvider.hpp>

namespace SpatialSEIR
{
    PerformDecorrelationStep::PerformDecorrelationStep(ModelContext* context_,
                                                       int iterationCount_)
    {
        context = new ModelContext*;
        iterationCount = new int;
        currentIteration = new int;

        *context = context_;
        *iterationCount = iterationCount_;
        *currentIteration = 0;
    }

    int PerformDecorrelationStep::getTaskType()
    {
        return(LSS_DECORRELATION_STEP_TASK_TYPE);
    }

    void PerformDecorrelationStep::executeTask()
    {
        int decorrIterationStride = ((*context) -> config -> useDecorrelation);
        if ((*currentIteration) % decorrIterationStride != 0)
        {
            *currentIteration += 1;
            return;
        }
        else
        {
            *currentIteration = 1;
        }


        int currentSamplingMode = (*((*context) -> beta_fc -> currentSampler)) -> getSamplerType(); 

        (*context) -> beta_fc -> evalCPU();
        (*context) -> beta_fc -> setSamplerType(PARAMETER_DECORR_SAMPLER);
        (*context) -> beta_fc -> evalCPU();



        int betaAccepted = *((*context) -> beta_fc -> accepted);
        int betaPrsAccepted = *((*context) -> betaPrs_fc -> accepted);

        *iterationCount = 0;
        double startingSliceWidth = *((*context) -> beta_fc -> sliceWidth);
        while (*((*context) -> beta_fc -> accepted) == betaAccepted &&
                *iterationCount < 1000)
        {
            (*context) -> beta_fc -> sample(0);
            (*((*context) -> beta_fc -> sliceWidth)) *= 0.99;
            *iterationCount += 1;
        }
        (*((*context) -> beta_fc -> sliceWidth)) = startingSliceWidth;

        if (*((*context) -> beta_fc -> accepted) == betaAccepted)
        {
            lssCout << "Decorrelation sampler did not update (beta).\n";
        }
        (*context) -> beta_fc -> setSamplerType(currentSamplingMode);

        (*context) -> betaPrs_fc -> setSamplerType(PARAMETER_DECORR_SAMPLER);
        if ((*context) -> config -> reinfectionMode != 1)
        {
            return;
        }

        *iterationCount = 0;
        while (*((*context) -> betaPrs_fc -> accepted) == betaPrsAccepted &&
                *iterationCount < 1000)
        {
            (*context) -> betaPrs_fc -> sample(0);
            *iterationCount += 1;
        }
        if (*((*context) -> betaPrs_fc -> accepted) == betaPrsAccepted)
        {
            lssCout << "Decorrelation sampler did not update (beta_pRS).\n";
        }

        (*context) -> betaPrs_fc -> setSamplerType(currentSamplingMode);
    }

    PerformDecorrelationStep::~PerformDecorrelationStep()
    {
        delete context;
        delete iterationCount;
        delete currentIteration;
    }
}

