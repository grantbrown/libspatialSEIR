#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<cstring>
#include<vector>
#endif

#ifndef LSS_SAMPLERS_INC
#define LSS_SAMPLERS_INC

#define COMPARTMENT_METROPOLIS_SAMPLER 1
#define COMPARTMENT_IDX_METROPOLIS_SAMPLER 2
#define COMPARTMENT_IDX_SLICE_SAMPLER 3
#define INITCOMPARTMENT_METROPOLIS_SAMPLER 4
#define INITCOMPARTMENT_IDX_METROPOLIS_SAMPLER 5
#define INITCOMPARTMENT_IDX_SLICE_SAMPLER 6
#define PARAMETER_SINGLE_METROPOLIS_SAMPLER 7
#define PARAMETER_JOINT_METROPOLIS_SAMPLER 8
#define PARAMETER_JOINT_SLICE_SAMPLER 9
#define PARAMETER_DECORR_SAMPLER 10

#define COMPARTMENT_METROPOLIS_SAMPLER_OCL 11
#define INITCOMPARTMENT_METROPOLIS_SAMPLER_OCL 12
#define PARAMETER_JOINT_METROPOLIS_SAMPLER_OCL 13

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
            virtual int getSamplerType() = 0;
    };

    // Samplers for CompartmentalModelMatrix classes
    class CompartmentMetropolisSampler : public Sampler
    {
        public: 
            CompartmentMetropolisSampler(ModelContext* context,
                                         CompartmentFullConditional* compartmentFC,
                                         int* compartmentData);
            void drawSample();
            int getSamplerType();
            ~CompartmentMetropolisSampler();

            ModelContext** context;
            CompartmentFullConditional** compartmentFC;
            int** compartmentData;
    };

    class CompartmentMetropolisSampler_OCL : public Sampler
    {
        public: 
            CompartmentMetropolisSampler_OCL(ModelContext* context,
                                         CompartmentFullConditional* compartmentFC,
                                         int* compartmentData);
            void drawSample();
            int getSamplerType();
            ~CompartmentMetropolisSampler_OCL();

            ModelContext** context;
            CompartmentFullConditional** compartmentFC;
            int** compartmentData;
    };

    class IndexedCompartmentMetropolisSampler : public Sampler
    {
        public: 
            IndexedCompartmentMetropolisSampler(ModelContext* context,
                                         CompartmentFullConditional* compartmentFC,
                                         int* compartmentData);
            void drawSample();
            int getSamplerType();
            ~IndexedCompartmentMetropolisSampler();

            ModelContext** context;
            CompartmentFullConditional** compartmentFC;
            int** indexLength;
            int** compartmentData;
            int** indexList;
    };

    class IndexedCompartmentSliceSampler : public Sampler
    {
        public: 
            IndexedCompartmentSliceSampler(ModelContext* context,
                                         CompartmentFullConditional* compartmentFC,
                                         int* compartmentData);
            void drawSample();
            int getSamplerType();
            ~IndexedCompartmentSliceSampler();

            ModelContext** context;
            CompartmentFullConditional** compartmentFC;
            int** indexLength;
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
            int getSamplerType();
            ~InitCompartmentMetropolisSampler();

            ModelContext** context;
            InitCompartmentFullConditional** initCompartmentFC;
            int** initCompartmentData;
    };

    class InitCompartmentMetropolisSampler_OCL : public Sampler
    {
        public:
            InitCompartmentMetropolisSampler_OCL(ModelContext* context,
                                             InitCompartmentFullConditional* initCompartmentFC,
                                             int* initCompartmentData);
            void drawSample();
            int getSamplerType();
            ~InitCompartmentMetropolisSampler_OCL();

            ModelContext** context;
            InitCompartmentFullConditional** initCompartmentFC;
            int** initCompartmentData;
    };

    class IndexedInitCompartmentMetropolisSampler : public Sampler
    {
        public:
            IndexedInitCompartmentMetropolisSampler(ModelContext* context,
                                             InitCompartmentFullConditional* initCompartmentFC,
                                             int* initCompartmentData);
            void drawSample();
            int getSamplerType();
            ~IndexedInitCompartmentMetropolisSampler();

            ModelContext** context;
            InitCompartmentFullConditional** initCompartmentFC;
            int** indexLength;
            int** indexList;
            int** initCompartmentData;

    };
    class IndexedInitCompartmentSliceSampler : public Sampler
    {
        public:
            IndexedInitCompartmentSliceSampler(ModelContext* context,
                                             InitCompartmentFullConditional* initCompartmentFC,
                                             int* initCompartmentData);
            void drawSample();
            int getSamplerType();
            ~IndexedInitCompartmentSliceSampler();

            ModelContext** context;
            InitCompartmentFullConditional** initCompartmentFC;
            int** indexLength;
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
            int getSamplerType();
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
            int getSamplerType();
            ~ParameterJointMetropolisSampler();

            double** param;
            ModelContext** context;
            ParameterFullConditional** paramFC;
    };

    class ParameterJointMetropolisSampler_OCL : public Sampler
    {
        public:
            ParameterJointMetropolisSampler_OCL(ModelContext* context,
                                            ParameterFullConditional* paramFC,
                                            double* param);
            void drawSample();
            int getSamplerType();
            ~ParameterJointMetropolisSampler_OCL();

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
            int getSamplerType();
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
            int getSamplerType();
            ~ParameterDecorrelationSampler();

            double** param;
            ModelContext** context;
            ParameterFullConditional** paramFC;
            CovariateMatrix** proposalMatrix;
    };

    









}
#endif
