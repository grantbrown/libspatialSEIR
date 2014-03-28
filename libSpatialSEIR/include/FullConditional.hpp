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

    class FullConditional
    {
        public:
            //Template for shared methods
            virtual ~FullConditional(){}; 
            virtual int cacheEvalCalculation(double* cachedValues) = 0;
            virtual int evalCPU() = 0;
            virtual int evalCPU(int startLoc, int startTime, double* cachedValues) = 0;
            virtual int evalOCL() = 0;
            virtual int sampleCPU() = 0;
            virtual int sampleOCL() = 0;
            virtual double getValue() = 0;
            // The following methods only apply to the compartment full conditionals, 
            // but require dummy implementations for all FullConditionals to keep 
            // the R dynamic loader happy. Consider factoring out separate FullConditional
            // superclasses in the future. 
            virtual int calculateRelevantCompartments() = 0;
            virtual int calculateRelevantCompartments(int startLoc, int startTime) = 0;

            //Declaration for inherited methods
            int sampleCompartment(ModelContext* context,
                                  InitData* A0,
                                  CompartmentalModelMatrix* drawCompartment,
                                  CompartmentalModelMatrix* destCompartment,
                                  CompartmentalModelMatrix* starCompartment,
                                  double width); 
    };

    class FC_S_Star : public FullConditional
    {
        public:
            FC_S_Star(ModelContext * _context,
                      CompartmentalModelMatrix *_S_star, 
                      CompartmentalModelMatrix *_S, 
                      CompartmentalModelMatrix *_R,
                      InitData *_A0,
                      CovariateMatrix *_X,
                      double *_p_se,
                      double *_p_rs,
                      double *_beta,
                      double *_rho);
            virtual int cacheEvalCalculation(double* cachedValues);
            virtual int evalCPU();
            virtual int evalCPU(int startLoc, int startTime, double* cachedValues);
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual double getValue();
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int startLoc, int startTime);
            virtual ~FC_S_Star();

            ModelContext **context;
            CompartmentalModelMatrix **S_star; 
            CompartmentalModelMatrix **S; 
            CompartmentalModelMatrix **R;
            InitData **A0;
            CovariateMatrix **X;
            double **p_se;
            double **p_rs;
            double **beta; 
            double **rho;
            double* value;

    };

    class FC_E_Star : public FullConditional
    {
        public:
            FC_E_Star(ModelContext *_context,
                      CompartmentalModelMatrix *_E_star, 
                      CompartmentalModelMatrix *_E, 
                      CompartmentalModelMatrix *_S, 
                      CovariateMatrix *_X,
                      InitData *_A0,
                      double *_p_se,
                      double *_p_ei,
                      double *_rho,
                      double *_beta);
            ~FC_E_Star();

            virtual int cacheEvalCalculation(double* cachedValues);
            virtual int evalCPU();
            virtual int evalCPU(int startLoc, int startTime, double* cachedValues);
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual double getValue();
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int startLoc, int startTime);

            ModelContext **context;
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **E; 
            CompartmentalModelMatrix **S; 
            CovariateMatrix **X;
            InitData **A0;
            double **p_se;
            double **p_ei;
            double **rho;
            double **beta;
            double* value;

    };

    class FC_R_Star : public FullConditional
    {
        public:
            FC_R_Star(ModelContext *_context,
                      CompartmentalModelMatrix *_R_star,
                      CompartmentalModelMatrix *_R,
                      CompartmentalModelMatrix *_I,
                      InitData *_A0,
                      double *_p_rs,
                      double *_p_ir);
            ~FC_R_Star();

            virtual int cacheEvalCalculation(double* cachedValues);
            virtual int evalCPU();
            virtual int evalCPU(int startLoc, int startTime, double* cachedValues);
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual double getValue();
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int startLoc, int startTime);

            ModelContext **context;
            CompartmentalModelMatrix **R_star;
            CompartmentalModelMatrix **R;
            CompartmentalModelMatrix **I;
            InitData **A0;
            double **p_rs;
            double **p_ir;
            double* value;

    };

    class FC_Beta : public FullConditional
    {
        public:
            FC_Beta(ModelContext *_context,
                    CompartmentalModelMatrix *_E_star, 
                    CompartmentalModelMatrix *_S, 
                    InitData *_A0,
                    CovariateMatrix *_X,
                    double *_p_se, 
                    double *_beta,
                    double *_rho); 
            ~FC_Beta();

            virtual int cacheEvalCalculation(double* cachedValues);
            virtual int evalCPU();
            virtual int evalCPU(int startLoc, int startTime, double* cachedValues);
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual double getValue();
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int startLoc, int startTime);

            ModelContext **context;
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **S; 
            InitData **A0;
            CovariateMatrix **X;
            double **p_se;
            double **beta;
            double **rho;
            double* value;

    };

    class FC_P_RS : public FullConditional
    {
        public:
            FC_P_RS(ModelContext *_context,
                    CompartmentalModelMatrix *_S_star,
                    CompartmentalModelMatrix *_R,
                    InitData *_A0,
                    double *_p_rs 
                    );
            ~FC_P_RS();
            virtual int cacheEvalCalculation(double* cachedValues);
            virtual int evalCPU();
            virtual int evalCPU(int startLoc, int startTime, double* cachedValues);
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual double getValue();
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int startLoc, int startTime);

            ModelContext **context;
            CompartmentalModelMatrix **S_star;
            CompartmentalModelMatrix **R;
            InitData **A0;
            double **p_rs;
            double* value;

    };

    class FC_Rho : public FullConditional 
    {
        public:
            FC_Rho(ModelContext *_context,
                   CompartmentalModelMatrix *_E_star, 
                   CompartmentalModelMatrix *_S, 
                   InitData *_A0,
                   CovariateMatrix *_X,
                   double *_p_se, 
                   double *_beta, 
                   double *_rho
                   );
            ~FC_Rho();
            virtual int cacheEvalCalculation(double* cachedValues);
            virtual int evalCPU();
            virtual int evalCPU(int startLoc, int startTime, double* cachedValues);
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual double getValue();
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int startLoc, int startTime);

            ModelContext **context;
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **S; 
            InitData **A0;
            CovariateMatrix **X;
            double **p_se;
            double **beta;
            double **rho;
            double* value;

    };
    

    class FC_P_EI : public FullConditional
    {
        public:
            FC_P_EI(ModelContext *_context,
                    CompartmentalModelMatrix *_I_star,
                    CompartmentalModelMatrix *_E,
                    InitData *_A0,
                    double *_p_ei);
            virtual int cacheEvalCalculation(double* cachedValues);
            virtual int evalCPU();
            virtual int evalCPU(int startLoc, int startTime, double* cachedValues);
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual double getValue();
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int startLoc, int startTime);

            ~FC_P_EI();
            ModelContext** context;
            CompartmentalModelMatrix **I_star;
            CompartmentalModelMatrix **E;
            InitData **A0;
            double **p_ei;
            double* value;
    };

    class FC_P_IR : public FullConditional
    {
        
        public:
            FC_P_IR(ModelContext *_context,
                    CompartmentalModelMatrix *_R_star,
                    CompartmentalModelMatrix *_I, 
                    InitData *_A0,
                    double *_p_ir);
            ~FC_P_IR();
            virtual int cacheEvalCalculation(double* cachedValues);
            virtual int evalCPU();
            virtual int evalCPU(int startLoc, int startTime, double* cachedValues);
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual double getValue();
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int startLoc, int startTime);

            ModelContext **context;
            CompartmentalModelMatrix **R_star;
            CompartmentalModelMatrix **I;
            InitData **A0;
            double **p_ir;
            double* value;
    };
}

#endif
