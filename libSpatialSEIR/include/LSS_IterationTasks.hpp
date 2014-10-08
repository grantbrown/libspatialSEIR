#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<cstring>
#include<vector>
#endif

#ifndef LSS_ITERATION_TASKS_INC
#define LSS_ITERATION_TASKS_INC

#define LSS_SAMPLING_INDEX_TASK_TYPE 1
#define LSS_DECORRELATION_STEP_TASK_TYPE 2
#define LSS_HYBRID_SE_EI_UPDATE_STEP_TASK_TYPE 3


namespace SpatialSEIR
{

    class ModelContext;
    class FC_Gamma_EI;
    class FC_Beta;
    class ParameterHybridSampler;

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
    class SetCompartmentSamplingIndicesTask : public IterationTask
    {
        public:
            SetCompartmentSamplingIndicesTask(ModelContext* context);
            ~SetCompartmentSamplingIndicesTask();
            void executeTask();
            int getTaskType();

            ModelContext** context;
            int** index;
            int** indexLength;
    };

    class PerformDecorrelationStep : public IterationTask
    {
        public:
            PerformDecorrelationStep(ModelContext* context, 
                                        int iterationCount);
            ~PerformDecorrelationStep();
            void executeTask();
            int getTaskType();

            ModelContext** context;
            int* iterationCount;
            int* currentIteration;
    };

    class PerformHybridSE_EI_UpdateStep : public IterationTask
    {
        public:
            PerformHybridSE_EI_UpdateStep(ModelContext* context,
                                          FC_Gamma_EI* fc_gammaEI,
                                          FC_Beta* fc_beta,
                                          int iterationCount);
            ~PerformHybridSE_EI_UpdateStep();
            void executeTask();
            int getTaskType();
            ModelContext** context;
            ParameterHybridSampler* sampler;
            int* iterationCount;
            int* currentIteration;
    };
}
#endif 
