#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<cstring>
#include<vector>
#endif

#include "OCLProvider.hpp"
#include "CompartmentalModelMatrix.hpp"
#include "CovariateMatrix.hpp"

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class InitData
    {
        public:
            InitData(double *_S0,
                        double *_E0,
                        double *_I0,
                        double *_R0,
                        double *_S_star0,
                        double *_E_star0,
                        double *_I_star0,
                        double *_R_star0);
            double *S0;
            double *E0;
            double *I0;
            double *R0;
            double *S_star0;
            double *E_star0;
            double *I_star0;
            double *R_star0;
    };

    class FullConditional
    {
        public:
            //Methods
            FullConditional();            
            ~FullConditional(); 
            int evalCPU();
            int evalOCL();
            int sampleCPU();
            int sampleOCL();
            //Attributes
            double *value;
    };

    class FC_S_Star : public FullConditional 
    {
        public:
            FC_S_Star(CompartmentalModelMatrix *_S_star, 
                      CompartmentalModelMatrix *_E_star, 
                      CompartmentalModelMatrix *_R_star,
                      InitData *_A0,
                      CovariateMatrix *_X,
                      double *_p_se,
                      double *_p_rs,
                      double *_beta,
                      double *_rho);
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();

            CompartmentalModelMatrix *S_star; 
            CompartmentalModelMatrix *E_star; 
            CompartmentalModelMatrix *R_star;
            InitData *A0;
            CovariateMatrix *X;
            double *p_se;
            double *p_rs;
            double *beta; 
            double *rho;
    };

    class FC_E_Star : public FullConditional 
    {
        public:
            FC_E_Star(CompartmentalModelMatrix *_S_star, 
                      CompartmentalModelMatrix *_E_star, 
                      CompartmentalModelMatrix *_I_star, 
                      CovariateMatrix *_X,
                      InitData *_A0,
                      double *_p_se,
                      double *_p_rs,
                      double *_rho,
                      double *_beta);

            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();

            CompartmentalModelMatrix *S_star; 
            CompartmentalModelMatrix *E_star; 
            CompartmentalModelMatrix *I_star; 
            CovariateMatrix *X;
            InitData *A0;
            double *p_se;
            double *p_rs;
            double *rho;
            double *beta;
    };

    class FC_R_Star : public FullConditional
    {
        public:
            FC_R_Star(CompartmentalModelMatrix *_R_star,
                      CompartmentalModelMatrix *_S_star,
                      CompartmentalModelMatrix *_I_star,
                      InitData *_A0,
                      double *_p_rs,
                      double *_p_ir);

            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();

            CompartmentalModelMatrix *R_star;
            CompartmentalModelMatrix *S_star;
            CompartmentalModelMatrix *I_star;
            InitData *A0;
            double* p_rs;
            double *p_ir;
    };

    class FC_Beta : public FullConditional 
    {
        public:
            FC_Beta(CompartmentalModelMatrix *_E_star, 
                    CompartmentalModelMatrix *_S_star, 
                    InitData *_A0,
                    CovariateMatrix *_X,
                    double *_p_se, 
                    double *_beta,
                    double *_rho); 

            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();

            CompartmentalModelMatrix *E_star; 
            CompartmentalModelMatrix *S_star; 
            InitData *A0;
            CovariateMatrix *X;
            double *p_se;
            double *beta;
            double *rho;
    };

    class FC_P_RS : public FullConditional
    {
        public:
            FC_P_RS(CompartmentalModelMatrix *_S_star,
                    CompartmentalModelMatrix *_R_star,
                    InitData *_A0,
                    double *_p_rs 
                    );
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();

            CompartmentalModelMatrix *S_star;
            CompartmentalModelMatrix *R_star;
            InitData *A0;
            double *p_rs;

    };

    class FC_Rho : public FullConditional
    {
        public:
            FC_Rho(CompartmentalModelMatrix *_S_star, 
                   CompartmentalModelMatrix *_E_star, 
                   InitData *_A0,
                   CovariateMatrix *_X,
                   double *_p_se, 
                   double *_beta, 
                   double *_rho
                   );
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();

            CompartmentalModelMatrix *S_star; 
            CompartmentalModelMatrix *E_star; 
            InitData *A0;
            CovariateMatrix *X;
            double *p_se;
            double *beta;
            double *rho;

    };
    

    class FC_P_EI : public FullConditional
    {
        public:
            FC_P_EI(CompartmentalModelMatrix *_I_star,
                    CompartmentalModelMatrix *_E_star,
                    InitData *_A0,
                    double *_p_ei);
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();

            CompartmentalModelMatrix *I_star;
            CompartmentalModelMatrix *E_star;
            InitData *A0;
            double *p_ei;
    };

    class FC_P_IR : public FullConditional 
    {
        
        public:
            FC_P_IR(CompartmentalModelMatrix *_I_star, 
                    CompartmentalModelMatrix *_R_star,
                    InitData *_A0,
                    double *_p_ir);
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();

            CompartmentalModelMatrix *I_star;
            CompartmentalModelMatrix *R_star;
            InitData *A0;
            double *p_ir;
    };

}
