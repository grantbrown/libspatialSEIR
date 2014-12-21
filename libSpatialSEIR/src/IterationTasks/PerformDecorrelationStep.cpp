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
#include<LSS_FC_Beta_P_RS.hpp>
#include<LSS_Samplers.hpp>
#include<IOProvider.hpp>

namespace SpatialSEIR
{
    PerformDecorrelationStep::PerformDecorrelationStep(ModelContext* context_,
                                                       int iterationCount_,
                                                       ParameterFullConditional* intensityFC_)
    {
        context = new ModelContext*;
        iterationCount = new int;
        currentIteration = new int;
        intensityFC = new ParameterFullConditional*;
        *context = context_;
        *iterationCount = iterationCount_;
        *currentIteration = 0;
        *intensityFC = intensityFC_;
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


        int currentSamplingMode = (*((*intensityFC) -> currentSampler)) -> getSamplerType(); 

        (*intensityFC) -> setSamplerType(PARAMETER_DECORR_SAMPLER);
        (*context) -> betaPrs_fc -> setSamplerType(PARAMETER_DECORR_SAMPLER);

        int betaAccepted = *((*intensityFC) -> accepted);
        int betaPrsAccepted = *((*context) -> betaPrs_fc -> accepted);
        double startingSliceWidth = *((*intensityFC) -> sliceWidth);
        *iterationCount = 0;
        if (*((*intensityFC) -> varLen) > 1)
        {
            while (*((*intensityFC) -> accepted) == betaAccepted &&
                    *iterationCount < 1000)
            {
                (*intensityFC) -> sample(0);
                (*((*intensityFC) -> sliceWidth)) *= 0.99;
                *iterationCount += 1;
            }
            (*((*intensityFC) -> sliceWidth)) = startingSliceWidth;

            if (*((*intensityFC) -> accepted) == betaAccepted && *((*context) -> config -> verbose))
            {
                lssCout << "Decorrelation sampler did not update (beta).\n";
            }

        }

        *iterationCount = 0;
        if (*((*context) -> betaPrs_fc -> varLen) > 1 && (*context) -> config -> reinfectionMode == 1)
        {
            while (*((*context) -> betaPrs_fc -> accepted) == betaPrsAccepted &&
                    *iterationCount < 1000)
            {
                (*context) -> betaPrs_fc -> sample(0);
                *iterationCount += 1;
            }
            if (*((*context) -> betaPrs_fc -> accepted) == betaPrsAccepted && *((*context) -> config -> verbose))
            {
                lssCout << "Decorrelation sampler did not update (beta_pRS).\n";
            }
        }

        (*intensityFC) -> setSamplerType(currentSamplingMode);
        (*context) -> betaPrs_fc -> setSamplerType(currentSamplingMode);
    }

    PerformDecorrelationStep::~PerformDecorrelationStep()
    {
        delete context;
        delete iterationCount;
        delete currentIteration;
    }
}

