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

namespace SpatialSEIR
{
    SetCompartmentSamplingIndicesTask::SetCompartmentSamplingIndicesTask(ModelContext* context_)
    {
        context = new ModelContext*;
        index = new int*;
        indexLength = new int*;
        *context = context_;
        *index = ((*context) -> indexList);
        *indexLength = ((*context) -> indexLength);
    }

    int SetCompartmentSamplingIndicesTask::getTaskType()
    {
        return(LSS_SAMPLING_INDEX_TASK_TYPE);
    }

    void SetCompartmentSamplingIndicesTask::executeTask()
    {
        int i;
        int len = **indexLength;
        int maxIdx = *((*context) -> S -> nrow)*(*((*context) -> S -> ncol)) - 1;
        for (i = 0; i < len; i++)
        {
            (*index)[i] = (*context)->random->uniform_int(0,maxIdx);
        }
    }

    SetCompartmentSamplingIndicesTask::~SetCompartmentSamplingIndicesTask()
    {
        delete context;
        delete index;
        delete indexLength;
    }
}

