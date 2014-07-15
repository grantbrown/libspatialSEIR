#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<cstring>
#include<vector>
#endif

#ifndef FULL_CONDITIONAL_INC
#define FULL_CONDITIONAL_INC

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;
    class Sampler;

    //! struct containing hyperparameters for beta, betaP_RS, P_EI, and P_IR
    struct priorControl
    {
        double betaPriorPrecision;
        double betaPrsPriorPrecision;
        double P_EI_priorAlpha;
        double P_EI_priorBeta;
        double P_IR_priorAlpha;
        double P_IR_priorBeta;
    };

    //! struct containing initial slice sampling tuning parameters. 
    struct sliceParameters
    {
        double* S_starWidth;
        double* E_starWidth;
        double* R_starWidth;
        double* S0Width;
        double* I0Width;
        double* betaWidth;
        double* betaPrsWidth;
        double* rhoWidth;
        double* gammaWidth;
        double* gammaEiWidth;
        double* gammaIrWidth;
    };

    //! Wrapper for cblas::dgemm
    int matMult(double* output, double * A, double * B, int Arow, int Acol, 
            int Brow, int Bcol, bool TransA, bool TransB, int ldA, int ldB, int ldC);


    //! Simple class containing the starting compartment sizes. 
    class InitData
    {
        public:
            InitData(int *_S0,
                        int *_E0,
                        int *_I0,
                        int *_R0, 
                        int *nLoc);
            InitData();
            void populate(int *_S0,
                          int *_E0,
                          int *_I0,
                          int *_R0,
                          int *nLoc);

            ~InitData();
            int *S0;
            int *E0;
            int *I0;
            int *R0;
            int *numLocations;
    };



    /**
     * The FullConditional class serves as the grandparent class for the various
     * full conditional distributions. This structure is helpful, because we can 
     * then implement general methods which apply to all child classes of FullConditional. 
     *
     */
    class FullConditional
    {
        public:
            //Template for shared methods
            virtual ~FullConditional(){}; 
            virtual void sample(int verbose) = 0;
            virtual int sampleCPU() = 0;
            virtual int sampleOCL() = 0;
            virtual long double getValue() = 0;
            virtual void setValue(long double value) = 0;
            virtual int calculateRelevantCompartments() = 0;
            virtual int calculateRelevantCompartments_OCL() = 0;
            virtual void updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange) = 0;
            double acceptanceRatio();
            double* sliceWidth;
            std::vector<Sampler*>* samplers;
            Sampler** currentSampler;
            int* useOCL;
            int* samples;
            int* accepted;
    };

    /**
     * The CompartmentFullConditional class inherits the structure of 
     * FullConditional, and provides general sampling methods which 
     * apply to all compartment full conditional distributions. These are:
     * 1. FC_S_Star, the transition compartment from recovered to susceptible
     * 2. FC_E_Star, the transition compartment from susceptible to exposed
     * 3. FC_R_Star, the transition compartment from infectious to recovered 
     */
    class CompartmentFullConditional : public FullConditional 
    {
        public:
            //Template for shared methods
            virtual ~CompartmentFullConditional(){}; 
            virtual int evalCPU() = 0;
            virtual int evalCPU(int startLoc, int startTime) = 0;
            virtual int evalOCL() = 0;
            virtual void sample(int verbose) = 0;
            virtual int sampleCPU() = 0;
            virtual int sampleOCL() = 0;
            virtual long double getValue() = 0;
            virtual void setValue(long double value) = 0;
            virtual int calculateRelevantCompartments() = 0;
            virtual int calculateRelevantCompartments_OCL() = 0;
            virtual int calculateRelevantCompartments(int startLoc, int startTime) = 0;
            virtual void printDebugInfo(int loc, int tpt) = 0;
            void updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange);

            /** Standard Metropolis proposal, centered at current value. */
            void proposeUpdate(int* initCompartment,
                               int initCompartmentLength,
                               double width,
                               double* metropolisComponent1,
                               double* metropolisComponent2);

            /** Partial Metropolis proposal, for indices in indexList. */
            void proposeUpdate(int* initCompartment,
                               int* indexList,
                               int indexListLength,
                               double width,
                               double* metropolisComponent1,
                               double* metropolisComponent2);

            int sampleCompartment_CPU(ModelContext* context,
                                  CompartmentalModelMatrix* destCompartment,
                                  double width); 

            int sampleEntireCompartment_CPU(ModelContext* context,
                                  CompartmentalModelMatrix* destCompartment,
                                  double width); 

            int sampleEntireCompartment2_CPU(ModelContext* context,
                                  CompartmentalModelMatrix* destCompartment,
                                  double width); 

            int sampleCompartment_OCL(ModelContext* context,
                                  CompartmentalModelMatrix* destCompartment,
                                  double width); 

            int sampleCompartmentLocation(int loc, ModelContext* context,
                                  CompartmentalModelMatrix* destCompartment,
                                  double width); 

            int sliceSampleCompartmentLocation(int loc, ModelContext* context,
                                               CompartmentalModelMatrix* destCompartment,
                                               double width);

            double* steadyStateConstraintPrecision;
    };

    /**
     * The ParameterFullConditional class inherits the structure of FullConditional,
     * and provides sampling methods for all of the non-compartment parameters. These
     * include the following:
     * 1. FC_Beta, the full conditional for the parameters controlling the exposure probability
     * 2. FC_Beta_P_RS, the full conditional for the parameters controlling the reinfection probability
     * 3. FC_Rho, the full conditional for the spatial dependence parameter
     * 4. FC_P_EI, the full conditional for the E to I transition probability
     * 5. FC_P_IR, the full conditional for the I to R transition probability 
     */
    class ParameterFullConditional : public FullConditional
    {
        public:
            //Template for shared methods
            virtual ~ParameterFullConditional(){}; 
            virtual int evalCPU() = 0;
            virtual int evalOCL() = 0;
            virtual void sample(int verbose) = 0;
            virtual int sampleCPU() = 0;
            virtual int sampleOCL() = 0;
            virtual long double getValue() = 0;
            virtual void setValue(long double value) = 0;
            virtual int calculateRelevantCompartments() = 0;
            virtual int calculateRelevantCompartments_OCL() = 0;
            virtual double acceptanceRatio(int i);
            void updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange);

            int *varLen;

            /** Standard Metropolis proposal, centered at current value. */
            void proposeUpdate(double* variable,
                                       int varLen,
                                       double* width,
                                       double* metropolisComponent1,
                                       double* metropolisComponent2);

            /** Decorrelation proposal */
            void proposeUpdate(double* variable,
                               int varLen,
                               CovariateMatrix* priorMatrix,
                               double* width,
                               double* metropolisComponent1,
                               double* metropolisComponent2);

            int sampleDouble(ModelContext* context, 
                             double* variable,
                             int varLen,
                             double* width);

            int sampleEntireDouble_CPU(ModelContext* context, 
                             double* variable,
                             int varLen,
                             double* width);

            int sampleEntireDouble_OCL(ModelContext* context, 
                             double* variable,
                             int varLen,
                             double* width);

            int sampleDouble_OCL(ModelContext* context, 
                             double* variable,
                             int varLen,
                             double* width);

            int sampleDoubleMetropolis(ModelContext* context, 
                                       double* variable,
                                       int varLen,
                                       double* width);
    };


    /**
     * The InitCompartmentFullConditional class inherits the structure of FullConditional,
     * and provides sampling methods for the initial compartment size parameters. These are:
     * 1. FC_S0, the initial susceptible group
     * 2. FC_E0, the initial exposed group
     * 3. FC_I0, the initial infectious group
     * 4. FC_R0, the initial recovered group
     */
    class InitCompartmentFullConditional : public FullConditional
    {
        public:

            virtual ~InitCompartmentFullConditional(){}; 
            virtual int evalCPU() = 0;
            virtual int evalCPU(int startLoc) = 0;
            virtual int evalOCL() = 0;
            virtual void sample(int verbose) = 0;
            virtual int sampleCPU() = 0;
            virtual int sampleOCL() = 0;
            virtual long double getValue() = 0;
            virtual void setValue(long double value) = 0;
            virtual int calculateRelevantCompartments() = 0;
            virtual int calculateRelevantCompartments_OCL() = 0;
            virtual int calculateRelevantCompartments(int startLoc) = 0;
            void updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange);
            virtual void printDebugInfo(int loc) = 0;

            /** Standard Metropolis proposal, centered at current value. */
            void proposeUpdate(int* initCompartment,
                               int initCompartmentLength,
                               double width,
                               double* metropolisComponent1,
                               double* metropolisComponent2);

            /** Partial Metropolis proposal, for indices in indexList. */
            void proposeUpdate(int* initCompartment,
                               int* indexList,
                               int indexListLength,
                               double width,
                               double* metropolisComponent1,
                               double* metropolisComponent2);

            int sampleCompartment_CPU(ModelContext* context,
                                      int* initCompartment, 
                                      double width); 

            int sampleEntireCompartment_CPU(ModelContext* context,
                                            int* initCompartment,
                                            double width); 

            int sampleEntireCompartment_OCL(ModelContext* context,
                                            int* initCompartment,
                                            double width); 

            int sampleCompartmentLocation(int loc, ModelContext* context,
                                          int* initCompartment,
                                          double width); 

    };


}

#endif
