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

    int matMult(double* output, double * A, double * B, int Arow, int Acol, int Brow, int Bcol, bool TransA, bool TransB );


    class InitData
    {
        public:
            InitData(int *_S0,
                        int *_E0,
                        int *_I0,
                        int *_R0,
                        int *_S_star0,
                        int *_E_star0,
                        int *_I_star0,
                        int *_R_star0,
                        int *nLoc);
            InitData();
            void populate(int *_S0,
                          int *_E0,
                          int *_I0,
                          int *_R0,
                          int *_S_star0,
                          int *_E_star0,
                          int *_I_star0,
                          int *_R_star0,
                          int *nLoc);

            ~InitData();
            int *S0;
            int *E0;
            int *I0;
            int *R0;
            int *S_star0;
            int *E_star0;
            int *I_star0;
            int *R_star0;
            int *numLocations;
    };

    // Full conditional distribution parent class
    class FullConditional
    {
        public:
            //Template for shared methods
            virtual ~FullConditional(){}; 
            virtual int evalCPU() = 0;
            virtual int evalOCL() = 0;
            virtual int sampleCPU() = 0;
            virtual int sampleOCL() = 0;
            virtual double getValue() = 0;
            virtual void setValue(double value) = 0;
            virtual int calculateRelevantCompartments() = 0;
            double* sliceWidth;
    };

    // Parent class for compartment full conditional distributions
    class CompartmentFullConditional : public FullConditional 
    {
        public:
            //Template for shared methods
            virtual ~CompartmentFullConditional(){}; 
            virtual int cacheEvalCalculation(double* cachedValues) = 0;
            virtual int cacheEvalCalculation(int* inStarCompartmentCache,
                                             int* outStarCompartmentCache,
                                             int* toCompartmentCache,
                                             int* fromCompartmentCache,
                                             double* likelihoodCache,
                                             double* steadyStateCache);
            virtual int updateEvalCache(int startLoc, int startTime, double* cachedValues) = 0;
            virtual int evalCPU() = 0;
            virtual int evalCPU(int startLoc, int startTime, double* cachedValues) = 0;
            virtual int evalOCL() = 0;
            virtual int sampleCPU() = 0;
            virtual int sampleOCL() = 0;
            virtual double getValue() = 0;
            virtual void setValue(double value) = 0;
            virtual int calculateRelevantCompartments() = 0;
            virtual int calculateRelevantCompartments(int startLoc, int startTime) = 0;

            //Declaration for inherited methods
            int sampleCompartment(ModelContext* context,
                                  InitData* A0, 
                                  CompartmentalModelMatrix* destCompartment,
                                  double width, double* compartmentCache); 



            void sampleCompartment2(ModelContext* context,
                                    InitData* A0,
                                    CompartmentalModelMatrix* inStarCompartment,
                                    CompartmentalModelMatrix* outStarCompartment,
                                    CompartmentalModelMatrix* toCompartment,
                                    CompartmentalModelMatrix* fromCompartment,
                                    double width,
                                    double* likelihoodCache,
                                    double* steadyStateCache,
                                    int* inStarCompartmentCache,
                                    int* outStarCompartmentCache,
                                    int* toCompartmentCache,
                                    int* fromCompartmentCache);

            void sampleCompartmentLocation(int* inStarVector,
                                           int* outStarVector,
                                           int* toCompartmentVector,
                                           int* fromCompartmentVector,
                                           double* likelihoodCacheVector,
                                           double* steadyStateCacheVector,
                                           int i,
                                           InitData* A0,
                                           double width,
                                           ModelContext* context);
            

            int sampleCompartmentMemoized(ModelContext* context,
                                            InitData* A0, 
                                            CompartmentalModelMatrix* destCompartment,
                                            int width, double* compartmentCache); 

            int sampleCompartmentMetropolis(ModelContext* context,
                                  InitData* A0, 
                                  CompartmentalModelMatrix* destCompartment,
                                  double width, double* compartmentCache); 

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
            virtual double getValue() = 0;
            virtual void setValue(double value) = 0;
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
            virtual int cacheEvalCalculation(double* cachedValues);
            virtual int updateEvalCache(int startLoc, int startTime, double* cachedValues);
            virtual int evalCPU();
            virtual int evalCPU(int startLoc, int startTime, double* cachedValues);
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual double getValue();
            virtual void setValue(double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int startLoc, int startTime);
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
            double* value;
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

            virtual int cacheEvalCalculation(double* cachedValues);
            virtual int updateEvalCache(int startLoc, int startTime, double* cachedValues);
            virtual int evalCPU();
            virtual int evalCPU(int startLoc, int startTime, double* cachedValues);
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual double getValue();
            virtual void setValue(double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int startLoc, int startTime);

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
            double* value;
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

            virtual int cacheEvalCalculation(double* cachedValues);
            virtual int updateEvalCache(int startLoc, int startTime, double* cachedValues);
            virtual int evalCPU();
            virtual int evalCPU(int startLoc, int startTime, double* cachedValues);
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual double getValue();
            virtual void setValue(double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int startLoc, int startTime);

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
            double* value;
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
            virtual double getValue();
            virtual void setValue(double val);
            virtual int calculateRelevantCompartments();

            ModelContext **context;
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **S; 
            InitData **A0;
            CovariateMatrix **X;
            double **p_se;
            double **beta;
            double **rho;
            double* value;
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
            virtual double getValue();
            virtual void setValue(double val);
            virtual int calculateRelevantCompartments();

            ModelContext **context;
            CompartmentalModelMatrix **S_star;
            CompartmentalModelMatrix **R;
            CovariateMatrix **X;
            InitData **A0;
            double **beta_p_rs;
            double **p_rs;
            double* tausq;
            double* value;
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
            virtual double getValue();
            virtual void setValue(double val);
            virtual int calculateRelevantCompartments();

            ModelContext **context;
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **S; 
            InitData **A0;
            CovariateMatrix **X;
            double **p_se;
            double **beta;
            double **rho;
            double* value;
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
            virtual double getValue();
            virtual void setValue(double val);
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
            double* value;
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
            virtual double getValue();
            virtual void setValue(double val);
            virtual int calculateRelevantCompartments();

            ~FC_P_EI();
            ModelContext** context;
            CompartmentalModelMatrix **I_star;
            CompartmentalModelMatrix **E;
            InitData **A0;
            double **p_ei;
            double* value;
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
            virtual double getValue();
            virtual void setValue(double val);
            virtual int calculateRelevantCompartments();

            ModelContext **context;
            CompartmentalModelMatrix **R_star;
            CompartmentalModelMatrix **I;
            InitData **A0;
            double **p_ir;
            double* value;
            double* priorAlpha;
            double* priorBeta;
            double* sliceWidth;
    };
}

#endif
