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

    struct gammaArgs
    {
        double* priorAlpha;
        double* priorBeta;
        double* gamma;
    };

    // Roll gamma prior stuff into
    // this struct?
    struct priorControl
    {
        double betaPriorPrecision;
        double betaPrsPriorPrecision;
        double P_EI_priorAlpha;
        double P_EI_priorBeta;
        double P_IR_priorAlpha;
        double P_IR_priorBeta;
    };

    struct sliceParameters
    {
        double* S_starWidth;
        double* E_starWidth;
        double* R_starWidth;
        double* betaWidth;
        double* betaPrsWidth;
        double* rhoWidth;
        double* gammaWidth;
    };

    int matMult(double* output, double * A, double * B, int Arow, int Acol, 
            int Brow, int Bcol, bool TransA, bool TransB, int ldA, int ldB, int ldC);


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

    // Full conditional distribution parent class
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
            double* sliceWidth;
    };

    // Parent class for compartment full conditional distributions
    class CompartmentFullConditional : public FullConditional 
    {
        public:
            //Template for shared methods
            virtual ~CompartmentFullConditional(){}; 
            virtual int evalCPU(int startLoc, int startTime) = 0;
            virtual int evalOCL() = 0;
            virtual int sampleCPU() = 0;
            virtual int sampleOCL() = 0;
            virtual long double getValue() = 0;
            virtual void setValue(long double value) = 0;
            virtual int calculateRelevantCompartments() = 0;
            virtual int calculateRelevantCompartments(int startLoc, int startTime) = 0;
            virtual void printDebugInfo(int loc, int tpt) = 0;

            int sampleCompartment_CPU(ModelContext* context,
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
            double* sliceWidth;
    };

    // Parent class for double precision scalar/vector full conditionals
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
            int sampleDouble(ModelContext* context, 
                             double* variable,
                             int varLen,
                             double width);
            int sampleDoubleMetropolis(ModelContext* context, 
                                       double* variable,
                                       int varLen,
                                       double width);

            double* sliceWidth;
    };


    // Parent class for S0-R0
    class InitCompartmentFullConditional : public FullConditional
    {
        public:

            virtual ~InitCompartmentFullConditional(){}; 
            virtual int evalCPU(int startLoc) = 0;
            virtual int evalOCL() = 0;
            virtual int sampleCPU() = 0;
            virtual int sampleOCL() = 0;
            virtual long double getValue() = 0;
            virtual void setValue(long double value) = 0;
            virtual int calculateRelevantCompartments() = 0;
            virtual int calculateRelevantCompartments(int startLoc) = 0;
            virtual void printDebugInfo(int loc) = 0;

            int sampleCompartment_CPU(ModelContext* context,
                                      int* initCompartment, 
                                      double width); 

            int sampleCompartmentLocation(int loc, ModelContext* context,
                                          int* initCompartment,
                                          double width); 
 
            double* sliceWidth;
    };

    class FC_S0 : InitCompartmentFullConditional
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
            virtual int evalCPU(int startLoc);
            virtual int evalOCL() ;
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double value);
            virtual int calculateRelevantCompartments();
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
            double* sliceWidth;
            long double* value;

    };


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
            virtual int evalCPU(int startLoc);
            virtual int evalOCL() ;
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double value);
            virtual int calculateRelevantCompartments();
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
            double* sliceWidth;
            long double* value;
    };

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
            virtual int evalCPU(int startLoc);
            virtual int evalOCL() ;
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double value);
            virtual int calculateRelevantCompartments();
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
            double* sliceWidth;
            long double* value;
    };

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
            virtual int evalCPU(int startLoc);
            virtual int evalOCL() ;
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double value);
            virtual int calculateRelevantCompartments();
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
            double* sliceWidth;
            long double* value;
    };




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
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
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
            double* sliceWidth;
    };

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
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
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
            double* sliceWidth;
    };

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
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
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
            double* sliceWidth;
    };

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

            ModelContext **context;
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **S; 
            InitData **A0;
            CovariateMatrix **X;
            double **p_se;
            double **beta;
            double **rho;
            long double* value;
            double* sliceWidth;
            double* priorPrecision;
    };

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

            ModelContext **context;
            CompartmentalModelMatrix **S_star;
            CompartmentalModelMatrix **R;
            CovariateMatrix **X;
            InitData **A0;
            double **beta_p_rs;
            double **p_rs;
            double* tausq;
            long double* value;
            double* sliceWidth;

    };

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

            ModelContext **context;
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **S; 
            InitData **A0;
            CovariateMatrix **X;
            double **p_se;
            double **beta;
            double **rho;
            long double* value;
            double* sliceWidth;

    };

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
            double* sliceWidth;

    };    

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

            ~FC_P_EI();
            ModelContext** context;
            CompartmentalModelMatrix **I_star;
            CompartmentalModelMatrix **E;
            InitData **A0;
            double **p_ei;
            long double* value;
            double* priorAlpha;
            double* priorBeta;
            double* sliceWidth;

    };

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

            ModelContext **context;
            CompartmentalModelMatrix **R_star;
            CompartmentalModelMatrix **I;
            InitData **A0;
            double **p_ir;
            long double* value;
            double* priorAlpha;
            double* priorBeta;
            double* sliceWidth;
    };
}

#endif
