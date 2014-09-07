#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<cstring>
#include<vector>
#endif

#ifndef LSS_SAMPLERS_INC
#define LSS_SAMPLERS_INC

#define COMPARTMENT_METROPOLIS_SAMPLER 1
#define COMPARTMENT_IDX_METROPOLIS_SAMPLER 2
#define COMPARTMENT_IDX_SLICE_SAMPLER 3
#define COMPARTMENT_BINOM_PROPOSAL_METROPOLIS_SAMPLER 14
#define COMPARTMENT_BINOM_PROPOSAL_SLICE_SAMPLER 15
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
    class FullConditional;
    class CompartmentFullConditional;
    class InitCompartmentFullConditional;
    class ParameterFullConditional;
    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class RandomNumberProvider;
    class OCLProvider;

    /** The sampler class is the parent class for all MCMC sampler objects, and guarantees that the same high level API
     * applies to all such objects.*/
    class Sampler
    {
        public:
            virtual ~Sampler(){};
            virtual void drawSample() = 0; 
            virtual int getSamplerType() = 0;
    };

    /** The CompartmentMetropolisSampler class is child of the Sampler class which draws samples from the 
     * posterior distribution of the various transition compartments using a full compartment Metropolis proposal. */
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

    /** The CompartmentMetropolisSampler_OCL class is child of the Sampler class which draws samples from the 
     * posterior distribution of the various transition compartments using a full compartment Metropolis proposal. This 
     * version of the Metropolis sampler uses the available OpenCL features of libspatialSEIR.*/
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

    /** The IndexedCompartmentMetropolisSampler class is child of the Sampler class which draws samples from the 
     * posterior distribution of the various transition compartments using a partial Metropolis proposal, with 
     * proposal indices determined by the ModelContext to which it belongs. This sampler is designed to reduce
     * autocorrelation by letting all compartments move together using the same indices. */
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
    /** The IndexedCompartmentSliceSampler class is child of the Sampler class which draws samples from the 
     * posterior distribution of the various transition compartments using a partial slice sampling proposal, with 
     * proposal indices determined by the ModelContext to which it belongs. This sampler is designed to reduce
     * autocorrelation by letting all compartments move together using the same indices. */
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


    /** The CompartmentBinomialMetropolisSampler class is child of the Sampler class which draws samples from the 
     * posterior distribution of the various transition compartments using a chain binomial proposal based on the parameters. */
    class CompartmentBinomialMetropolisSampler : public Sampler
    {
        public: 
            CompartmentBinomialMetropolisSampler(ModelContext* context,
                                         CompartmentFullConditional* compartmentFC,
                                         int* compartmentData,
                                         int* compartmentFrom, 
                                         int* compartmentTo,
                                         double* probabilityVector,
                                         int probabilityVectorLen);
            void drawSample();
            int getSamplerType();
            void genProposal();
            ~CompartmentBinomialMetropolisSampler();

            ModelContext** context;
            CompartmentFullConditional** compartmentFC;
            int** compartmentData;
            int** compartmentFrom;
            int** compartmentTo;
            double** probabilityVector;
            int* probabilityVectorLen;
    };


    /** The CompartmentBinomialSliceSampler class is child of the Sampler class which draws samples from the 
     * posterior distribution of the various transition compartments using a chain binomial proposal based on the parameters. */
    class CompartmentBinomialSliceSampler : public Sampler
    {
        public: 
            CompartmentBinomialSliceSampler(ModelContext* context,
                                         CompartmentFullConditional* compartmentFC,
                                         int* compartmentData,
                                         int* compartmentFrom, 
                                         int* compartmentTo,
                                         double* probabilityVector,
                                         int probabilityVectorLen);
            void drawSample();
            int getSamplerType();
            void genProposal();
            ~CompartmentBinomialSliceSampler();

            ModelContext** context;
            CompartmentFullConditional** compartmentFC;
            int** compartmentData;
            int** compartmentFrom;
            int** compartmentTo;
            double** probabilityVector;
            int* probabilityVectorLen;
    };





    /** The InitCompartmentMetropolisSampler functions identically to the CompartmentMetropolisSampler class, but for 
     * the vector of initial values of each CompartmentalModelMatrix.*/
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

    /** The InitCompartmentMetropolisSampler_OCL functions identically to the CompartmentMetropolisSampler_OCL class, but for 
     * the vector of initial values of each CompartmentalModelMatrix. This class makes use of the OpenCL functionality of libspatialSEIR*/
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

    /** The IndexedInitCompartmentMetropolisSampler is depricated.*/
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

    /** The IndexedInitCompartmentSliceSampler is depricated.*/
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
    
    /** The ParameterSingleMetropolisSampler uses Metropolis proposals of size one to 
     * sample non CompartmentalMatrix parameters. 
     * */
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

     /** The ParameterJointMetropolisSampler uses full vector Metropolis proposals to 
     * sample non CompartmentalMatrix parameters. 
     * */
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

     /** The ParameterJointMetropolisSampler_OCL uses full vector Metropolis proposals to 
     * sample non CompartmentalMatrix parameters. This class makes use of the OpenCL Capabilities
     * of libspatialSEIR.
     * */
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

     /** The ParameterJointSliceSampler uses full vector slice sampling proposals to 
     * sample non CompartmentalMatrix parameters. 
     * */
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

    /** The decorrelation sampler functions identically to the ParameterJointMetropolisSampler class, 
     * except that proposals are drawn from the null space of the relevant design matrix. This can
     * be useful when the explanatory variables are correlated.*/
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
            double* proposalCache;
            double* proposalCache2;
            ModelContext** context;
            ParameterFullConditional** paramFC;
            CovariateMatrix** proposalMatrix;
    }; 

}
#endif
