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

    //! DEPRICATED - Prior arguments for the gamma term (external infection probability) 
    struct gammaArgs
    {
        double* priorAlpha;
        double* priorBeta;
        double* gamma;
    };

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
            virtual int sampleCPU() = 0;
            virtual int sampleOCL() = 0;
            virtual long double getValue() = 0;
            virtual void setValue(long double value) = 0;
            virtual int calculateRelevantCompartments() = 0;
            virtual int calculateRelevantCompartments_OCL() = 0;
            void updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange);
            double acceptanceRatio();
            double* sliceWidth;
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
            virtual int sampleCPU() = 0;
            virtual int sampleOCL() = 0;
            virtual long double getValue() = 0;
            virtual void setValue(long double value) = 0;
            virtual int calculateRelevantCompartments() = 0;
            virtual int calculateRelevantCompartments_OCL() = 0;
            virtual int calculateRelevantCompartments(int startLoc, int startTime) = 0;
            virtual void printDebugInfo(int loc, int tpt) = 0;

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

            int jointSampleCompartmentLocation(int loc, int batchSize, int numBatches, int* batchCache, ModelContext* context,
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
            virtual int sampleCPU() = 0;
            virtual int sampleOCL() = 0;
            virtual long double getValue() = 0;
            virtual void setValue(long double value) = 0;
            virtual int calculateRelevantCompartments() = 0;
            virtual int calculateRelevantCompartments_OCL() = 0;
            int sampleDouble(ModelContext* context, 
                             double* variable,
                             int varLen,
                             double width);

            int sampleEntireDouble_CPU(ModelContext* context, 
                             double* variable,
                             int varLen,
                             double width);

            int sampleEntireDouble_OCL(ModelContext* context, 
                             double* variable,
                             int varLen,
                             double width);

            int sampleDouble_OCL(ModelContext* context, 
                             double* variable,
                             int varLen,
                             double width);

            int sampleDoubleMetropolis(ModelContext* context, 
                                       double* variable,
                                       int varLen,
                                       double width);
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
            virtual int sampleCPU() = 0;
            virtual int sampleOCL() = 0;
            virtual long double getValue() = 0;
            virtual void setValue(long double value) = 0;
            virtual int calculateRelevantCompartments() = 0;
            virtual int calculateRelevantCompartments_OCL() = 0;
            virtual int calculateRelevantCompartments(int startLoc) = 0;
            virtual void printDebugInfo(int loc) = 0;

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

    /**
     * FC_S0 gives the full conditional distribution for the vector of initially
     * susceptible individuals. 
     */
    class FC_S0 : public InitCompartmentFullConditional
    {
        public:

            FC_S0(ModelContext *_context,
                      CompartmentalModelMatrix *_S, 
                      CompartmentalModelMatrix *_E, 
                      CompartmentalModelMatrix *_E_star,
                      CompartmentalModelMatrix *_I_star, 
                      InitData *_A0,
                      double *_p_se,
                      double *_p_ei,
                      double sliceWidth);
            virtual ~FC_S0();
            virtual int evalCPU();
            virtual int evalCPU(int startLoc);
            virtual int evalOCL() ;
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double value);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();
            virtual int calculateRelevantCompartments(int startLoc);
            virtual void printDebugInfo(int loc);


            ModelContext** context;
            CompartmentalModelMatrix** S;
            CompartmentalModelMatrix** E;
            CompartmentalModelMatrix** E_star;
            CompartmentalModelMatrix** I_star;
            InitData** A0;
            double** p_se; 
            double** p_ei;
            long double* value;

    };

    /**
     * FC_E0 gives the full conditional distribution for the vector of initially
     * exposed individuals. 
     */
    class FC_E0 : public InitCompartmentFullConditional
    { 
        public:
            FC_E0(ModelContext *_context,
                      CompartmentalModelMatrix *_S, 
                      CompartmentalModelMatrix *_E, 
                      CompartmentalModelMatrix *_I, 
                      CompartmentalModelMatrix *_E_star,
                      CompartmentalModelMatrix *_I_star,
                      CompartmentalModelMatrix *_R_star, 
                      InitData *_A0,
                      double *_p_ir,
                      double *_p_ei,
                      double *_p_se,
                      double sliceWidth);
            virtual ~FC_E0(); 
            virtual int evalCPU();
            virtual int evalCPU(int startLoc);
            virtual int evalOCL() ;
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double value);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();
            virtual int calculateRelevantCompartments(int startLoc);
            virtual void printDebugInfo(int loc);

            ModelContext** context;
            CompartmentalModelMatrix** S;
            CompartmentalModelMatrix** E;
            CompartmentalModelMatrix** I;
            CompartmentalModelMatrix** E_star;
            CompartmentalModelMatrix** I_star;
            CompartmentalModelMatrix** R_star;
            InitData** A0;
            double** p_se;
            double** p_ir;
            double** p_ei;
            long double* value;
    };

    /**
     *
     * FC_I0 gives the full conditional distribution for the vector of initially
     * infectious individuals. 
     */
    class FC_I0 : public InitCompartmentFullConditional
    {
        public:
            FC_I0(ModelContext *_context,
                      CompartmentalModelMatrix *_S, 
                      CompartmentalModelMatrix *_I, 
                      CompartmentalModelMatrix *_R, 
                      CompartmentalModelMatrix *_S_star,
                      CompartmentalModelMatrix *_E_star,
                      CompartmentalModelMatrix *_R_star,
                      InitData *_A0,
                      double *_p_ir,
                      double *_p_rs,
                      double *_p_se,
                      double sliceWidth);
            virtual ~FC_I0(); 
            virtual int evalCPU();
            virtual int evalCPU(int startLoc);
            virtual int evalOCL() ;
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double value);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();
            virtual int calculateRelevantCompartments(int startLoc);
            virtual void printDebugInfo(int loc);

            ModelContext** context;
            CompartmentalModelMatrix** S;
            CompartmentalModelMatrix** I;
            CompartmentalModelMatrix** R;
            CompartmentalModelMatrix** S_star;
            CompartmentalModelMatrix** E_star;
            CompartmentalModelMatrix** R_star;
            InitData** A0;
            double** p_ir;
            double** p_rs;
            double** p_se; 
            long double* value;
    };

    /**
     * FC_R0 gives the full conditional distribution for the vector of initially
     * removed/recovered individuals. 
     */
    class FC_R0 : public InitCompartmentFullConditional
    {
        public:
            FC_R0(ModelContext *_context,
                      CompartmentalModelMatrix *_R, 
                      CompartmentalModelMatrix *_S,
                      CompartmentalModelMatrix *_S_star, 
                      CompartmentalModelMatrix *_E_star, 
                      CompartmentalModelMatrix *_R_star,
                      InitData *_A0,
                      double *_p_rs,
                      double *_p_se,
                      double sliceWidth);
            virtual ~FC_R0(); 
            virtual int evalCPU();
            virtual int evalCPU(int startLoc);
            virtual int evalOCL() ;
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double value);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();
            virtual int calculateRelevantCompartments(int startLoc);
            virtual void printDebugInfo(int loc);

            ModelContext** context;
            CompartmentalModelMatrix** R;
            CompartmentalModelMatrix** S;
            CompartmentalModelMatrix** S_star;
            CompartmentalModelMatrix** E_star; 
            CompartmentalModelMatrix** R_star;
            InitData** A0;
            double** p_rs;
            double** p_se; 
            long double* value;
    };




    /**
     * FC_S_Star gives the full conditional distribution of S_star, the 
     * transition matrix capturing individuals moving from the removed/recovered
     * category to the susceptible category. 
     */
    class FC_S_Star : public CompartmentFullConditional
    {
        public:
            FC_S_Star(ModelContext * _context,
                      CompartmentalModelMatrix *_S_star, 
                      CompartmentalModelMatrix *_S, 
                      CompartmentalModelMatrix *_R,
                      CompartmentalModelMatrix *_E_star,
                      CompartmentalModelMatrix *_R_star,
                      InitData *_A0,
                      CovariateMatrix *_X,
                      double *_p_se,
                      double *_p_rs,
                      double *_beta,
                      double *_rho,
                      double _steadyStateConstraintPrecision,
                      double sliceWidth);
            virtual int evalCPU(int startLoc, int startTime);
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();
            virtual int calculateRelevantCompartments(int startLoc, int startTime);
            virtual void printDebugInfo(int loc, int tpt);
            virtual ~FC_S_Star();

            ModelContext **context;
            CompartmentalModelMatrix **S_star; 
            CompartmentalModelMatrix **S; 
            CompartmentalModelMatrix **R;
            CompartmentalModelMatrix **E_star;
            CompartmentalModelMatrix **R_star;
            InitData **A0;
            CovariateMatrix **X;
            double **p_se;
            double **p_rs;
            double **beta; 
            double **rho;
            long double* value;
            double* steadyStateConstraintPrecision;
    };

    /**
     * FC_E_Star gives the full conditional distribution of E_star, the 
     * transition matrix capturing individuals moving from the susceptible
     * category to the exposed category. 
     */
    class FC_E_Star : public CompartmentFullConditional
    {
        public:
            FC_E_Star(ModelContext *_context,
                      CompartmentalModelMatrix *_E_star, 
                      CompartmentalModelMatrix *_E, 
                      CompartmentalModelMatrix *_S, 
                      CompartmentalModelMatrix *_I_star,
                      CovariateMatrix *_X,
                      InitData *_A0,
                      double *_p_se,
                      double *_p_ei,
                      double *_rho,
                      double *_beta,
                      double _steadyStateConstraintPrecision,
                      double sliceWidth);
            ~FC_E_Star();

            virtual int evalCPU(int startLoc, int startTime);
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();
            virtual int calculateRelevantCompartments(int startLoc, int startTime);
            virtual void printDebugInfo(int loc, int tpt);
            ModelContext **context;
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **E; 
            CompartmentalModelMatrix **S; 
            CompartmentalModelMatrix **I_star;
            CovariateMatrix **X;
            InitData **A0;
            double **p_se;
            double **p_ei;
            double **rho;
            double **beta;
            long double* value;
            double* steadyStateConstraintPrecision;
    };

    /**
     * FC_R_Star gives the full conditional distribution of R_star, the 
     * transition matrix capturing individuals moving from the infectious
     * category to the recovered/removed category. 
     */
    class FC_R_Star : public CompartmentFullConditional
    {
        public:
            FC_R_Star(ModelContext *_context,
                      CompartmentalModelMatrix *_R_star,
                      CompartmentalModelMatrix *_R,
                      CompartmentalModelMatrix *_I,
                      CompartmentalModelMatrix *_S_star,
                      CompartmentalModelMatrix *_E_star,
                      CompartmentalModelMatrix *_I_star,
                      CompartmentalModelMatrix *_S,
                      InitData *_A0,
                      double *_p_rs,
                      double *_p_ir,
                      double *_p_se,
                      double _steadyStateConstraintPrecision,
                      double sliceWidth);
            ~FC_R_Star();

            virtual int evalCPU(int startLoc, int startTime);
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();
            virtual int calculateRelevantCompartments(int startLoc, int startTime);
            virtual void printDebugInfo(int loc, int tpt);
            ModelContext **context;
            CompartmentalModelMatrix **R_star;
            CompartmentalModelMatrix **R;
            CompartmentalModelMatrix **I;
            CompartmentalModelMatrix **S_star;
            CompartmentalModelMatrix **E_star;
            CompartmentalModelMatrix **I_star;
            CompartmentalModelMatrix **S;
            InitData **A0;
            double **p_rs;
            double **p_ir;
            double **p_se;
            long double* value;
            double* steadyStateConstraintPrecision;
    };

    /**
     * FC_Beta gives the full conditional distribution of beta, the 
     * vector of regression parameters capturing the exposure intensity 
     * process. 
     */
    class FC_Beta : public ParameterFullConditional
    {
        public:
            FC_Beta(ModelContext *_context,
                    CompartmentalModelMatrix *_E_star, 
                    CompartmentalModelMatrix *_S, 
                    InitData *_A0,
                    CovariateMatrix *_X,
                    double *_p_se, 
                    double *_beta,
                    double *_rho,
                    double sliceWidth,
                    double _priorPrecision); 
            ~FC_Beta();

            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();

            ModelContext **context;
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **S; 
            InitData **A0;
            CovariateMatrix **X;
            double **p_se;
            double **beta;
            double **rho;
            long double* value;
            double* priorPrecision;
    };

    /**
     * FC_Beta_P_RS gives the full conditional distribution of beta_p_rs, the 
     * vector of regression parameters capturing the probability that an individual
     * transitions from R to S, the removed/recovered category to the susceptible category. 
     */
    class FC_Beta_P_RS : public ParameterFullConditional
    {
        public:
            FC_Beta_P_RS(ModelContext *_context,
                    CompartmentalModelMatrix *_S_star,
                    CompartmentalModelMatrix *_R,
                    CovariateMatrix* _X,
                    InitData *_A0,
                    double *_p_rs,
                    double *_beta_p_rs,
                    double _tausq,
                    double _sliceWidth
                    );
            ~FC_Beta_P_RS();
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();

            ModelContext **context;
            CompartmentalModelMatrix **S_star;
            CompartmentalModelMatrix **R;
            CovariateMatrix **X;
            InitData **A0;
            double **beta_p_rs;
            double **p_rs;
            double* tausq;
            long double* value;
    };

    /**
     * FC_Rho gives the full conditional distribution of rho, the 
     * scalar spatial dependence parameter. 
     */
    class FC_Rho : public ParameterFullConditional 
    {
        public:
            FC_Rho(ModelContext *_context,
                   CompartmentalModelMatrix *_E_star, 
                   CompartmentalModelMatrix *_S, 
                   InitData *_A0,
                   CovariateMatrix *_X,
                   double *_p_se, 
                   double *_beta, 
                   double *_rho,
                   double sliceWidth
                   );
            ~FC_Rho();
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();

            ModelContext **context;
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **S; 
            InitData **A0;
            CovariateMatrix **X;
            double **p_se;
            double **beta;
            double **rho;
            long double* value;
    };

    /**
     * FC_Gamma is depricated. Gamma was used to describe a time varying external 
     * infection source, but was found to be unneccesary for most epidemic models, 
     * as well as a potential source for identifiability issues. The code remains 
     * in case this decision changes. 
     */
    class FC_Gamma : public ParameterFullConditional 
    {
        public:
            FC_Gamma(ModelContext *_context,
                   CompartmentalModelMatrix *_E_star, 
                   CompartmentalModelMatrix *_S, 
                   InitData *_A0,
                   CovariateMatrix *_X,
                   double *_p_se, 
                   double *_beta, 
                   double *_gamma,
                   double *_priorAlpha,
                   double *_priorBeta,
                   double sliceWidth
                   );
            ~FC_Gamma();
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();

            ModelContext **context;
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **S; 
            InitData **A0;
            CovariateMatrix **X;
            double **p_se;
            double **beta;
            double **gamma;
            double* priorAlpha;
            double* priorBeta;
            long double* value;
    };    

    /**
     * FC_P_EI gives the full conditional distribution of p_ei, the 
     * probability that an exposed individual becomes infectious at a
     * given time point. 
     */
    class FC_P_EI : public ParameterFullConditional
    {
        public:
            FC_P_EI(ModelContext *_context,
                    CompartmentalModelMatrix *_I_star,
                    CompartmentalModelMatrix *_E,
                    InitData *_A0,
                    double *_p_ei,
                    double _priorAlpha,
                    double _priorBeta);
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();

            ~FC_P_EI();
            ModelContext** context;
            CompartmentalModelMatrix **I_star;
            CompartmentalModelMatrix **E;
            InitData **A0;
            double **p_ei;
            long double* value;
            double* priorAlpha;
            double* priorBeta;
    };

    /**
     * FC_P_IR gives the full conditional distribution of p_ir, the 
     * probability that an infectious individual recovers/is removed at 
     * a given time point.
     */
    class FC_P_IR : public ParameterFullConditional
    {
        
        public:
            FC_P_IR(ModelContext *_context,
                    CompartmentalModelMatrix *_R_star,
                    CompartmentalModelMatrix *_I, 
                    InitData *_A0,
                    double *_p_ir,
                    double _priorAlpha,
                    double _priorBeta);
            ~FC_P_IR();
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();

            ModelContext **context;
            CompartmentalModelMatrix **R_star;
            CompartmentalModelMatrix **I;
            InitData **A0;
            double **p_ir;
            long double* value;
            double* priorAlpha;
            double* priorBeta;
    };
}

#endif
