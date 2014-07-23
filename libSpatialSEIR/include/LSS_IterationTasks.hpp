#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<cstring>
#include<vector>
#endif

#ifndef LSS_ITERATION_TASKS_INC
#define LSS_ITERATION_TASKS_INC

#define LSS_SAMPLING_INDEX_TASK_TYPE 1

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class ModelContext;

    /** IterationTask instances are bits of logical code that need to 
     * execute periodically during MCMC sampling and do not fall into 
     * the usual FullConditional workflow. */
    class IterationTask
    {
        public:
            virtual ~IterationTask(){};
            virtual void executeTask() = 0; 
            virtual int getTaskType() = 0;
    };

    /** Several of the compartment sampling classes rely on a list 
     * of indices managed by model context to jointly update the 
     * same space-time locations for each of the compartments in order
     * to reduce sampler autocorrelation. This task updates the list of 
     * indices after each MCMC sample. */
    class SetCompartmentSamplingIndices : public IterationTask
    {
        public:
            SetCompartmentMetropolisIndices(ModelContext* context, 
                                            int updateProportion);
            void executeTask();
            void getTaskType();
    };
}
 
