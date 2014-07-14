#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<cstring>
#include<vector>
#endif

#ifndef LSS_SAMPLERS_INC
#define LSS_SAMPLERS_INC

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class FullConditional;
    class CompartmentFullConditional;
    class InitCompartmentFullConditional;
    class ParameterFullConditional;
    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class RandomNumberProvider;
    class OCLProvider;

    class Sampler
    {
        public:
            virtual ~Sampler(){};
            virtual void drawSample() = 0; 
    };

    // Samplers for CompartmentalModelMatrix classes
    class CompartmentMetropolisSampler : public Sampler
    {
        public: 
            CompartmentMetropolisSampler(ModelContext* context,
                                         CompartmentFullConditional* compartmentFC,
                                         int* compartmentData);
            void drawSample();
            ~CompartmentMetropolisSampler();

            ModelContext** context;
            CompartmentFullConditional** compartmentFC;
            int** compartmentData;
    };

    class IndexedCompartmentMetropolisSampler : public Sampler
    {
        public: 
            IndexedCompartmentMetropolisSampler(ModelContext* context,
                                         CompartmentFullConditional* compartmentFC,
                                         int* compartmentData,
                                         int* indexList,
                                         int indexLength);
            void drawSample();
            ~IndexedCompartmentMetropolisSampler();

            ModelContext** context;
            CompartmentFullConditional** compartmentFC;
            int* indexLength;
            int** compartmentData;
            int** indexList;
    };

    class IndexedCompartmentSliceSampler : public Sampler
    {
        public: 
            IndexedCompartmentSliceSampler(ModelContext* context,
                                         CompartmentFullConditional* compartmentFC,
                                         int* compartmentData,
                                         int* indexList,
                                         int indexLength);
            void drawSample();
            ~IndexedCompartmentSliceSampler();

            ModelContext** context;
            CompartmentFullConditional** compartmentFC;
            int* indexLength;
            int** compartmentData;
            int** indexList;
    };

    // Samplers for InitCompartment classes 
    class InitCompartmentMetropolisSampler : public Sampler
    {
        public:
            InitCompartmentMetropolisSampler(ModelContext* context,
                                             InitCompartmentFullConditional* initCompartmentFC,
                                             int* initCompartmentData);
            void drawSample();
            ~InitCompartmentMetropolisSampler();

            ModelContext** context;
            InitCompartmentFullConditional** initCompartmentFC;
            int** initCompartmentData;
    };
    class IndexedInitCompartmentMetropolisSampler : public Sampler
    {
        public:
            IndexedInitCompartmentMetropolisSampler(ModelContext* context,
                                             InitCompartmentFullConditional* initCompartmentFC,
                                             int* initCompartmentData,
                                             int* indexList,
                                             int indexLength);
            void drawSample();
            ~IndexedInitCompartmentMetropolisSampler();

            ModelContext** context;
            InitCompartmentFullConditional** initCompartmentFC;
            int* indexLength;
            int** indexList;
            int** initCompartmentData;

    };
    class IndexedInitCompartmentSliceSampler : public Sampler
    {
        public:
            IndexedInitCompartmentSliceSampler(ModelContext* context,
                                             InitCompartmentFullConditional* initCompartmentFC,
                                             int* initCompartmentData,
                                             int* indexList,
                                             int indexLength);
            void drawSample();
            ~IndexedInitCompartmentSliceSampler();

            ModelContext** context;
            InitCompartmentFullConditional** initCompartmentFC;
            int* indexLength;
            int** indexList;
            int** initCompartmentData;
    };
    
    // Samplers for non-compartmental parameters.
    class ParameterSingleMetropolisSampler : public Sampler
    {
        public:
            ParameterSingleMetropolisSampler(ModelContext* context,
                                            ParameterFullConditional* paramFC,
                                            double* param);
            void drawSample();
            ~ParameterSingleMetropolisSampler();

            double** param;
            ModelContext** context;
            ParameterFullConditional** paramFC;
    };
    
    class ParameterJointMetropolisSampler : public Sampler
    {
        public:
            ParameterJointMetropolisSampler(ModelContext* context,
                                            ParameterFullConditional* paramFC,
                                            double* param);
            void drawSample();
            ~ParameterJointMetropolisSampler();

            double** param;
            ModelContext** context;
            ParameterFullConditional** paramFC;
    };

    class ParameterJointSliceSampler : public Sampler
    {
        public:
            ParameterJointSliceSampler(ModelContext* context,
                                            ParameterFullConditional* paramFC,
                                            double* param);
            void drawSample();
            ~ParameterJointSliceSampler();

            double** param;
            ModelContext** context;
            ParameterFullConditional** paramFC;
    };

    class ParameterDecorrelationSampler : public Sampler
    {
        public:
            ParameterDecorrelationSampler(ModelContext* context,
                                          ParameterFullConditional* paramFC,
                                          double* param,
                                          CovariateMatrix* proposalMatrix);
            void drawSample();
            ~ParameterDecorrelationSampler();

            double** param;
            ModelContext** context;
            ParameterFullConditional** paramFC;
            CovariateMatrix** proposalMatrix;
    };

    









}
#endif
