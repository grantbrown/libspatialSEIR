#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<cstring>
#include<vector>
#endif

#include "OCLProvider.hpp"
#include "DataStructures/CompartmentalModelMatrix.hpp"
#include "DataStructures/CovariateMatrix.hpp"

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class InitData
    {
        public:
            InitData(double *S0,
                        double *E0,
                        double *I0,
                        double *R0,
                        double *S_star0,
                        double *E_star0,
                        double *I_star0,
                        double *R_star0);
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
            FC_S_Star(CompartmentalModelMatrix *S_star, 
                      CompartmentalModelMatrix *E_star, 
                      CompartmentalModelMatrix *R_star,
                      InitData *A0,
                      CovariateMatrix *X,
                      double *p_se,
                      double *p_rs,
                      double *beta,
                      double *rho);

            CompartmentalModelMatrix *S_star; 
            CompartmentalModelMatrix *E_star; 
            CompartmentalModelMatrix *R_star;
            InitData *A0;
            CovariateMatrix *X;
            double *p_se;
            double *p_rs,
            double *beta, 
            double *rho;
    };

    class FC_E_Star : public FullConditional 
    {
        public:
            FC_E_Star(CompartmentalModelMatrix *S_star, 
                      CompartmentalModelMatrix *E_star, 
                      CompartmentalModelMatrix *I_star, 
                      CovariateMatrix *X,
                      InitData *A0,
                      double *p_se,
                      double *p_rs,
                      double *rho,
                      double *beta);

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
            FC_R_Star(CompartmentalModelMatrix *R_Star,
                      CompartmentalModelMatrix *S_Star,
                      CompartmentalModelMatrix *I_star,
                      InitData *A0,
                      double* p_rs,
                      double *p_ir);

            CompartmentalModelMatrix *R_Star;
            CompartmentalModelMatrix *S_Star;
            CompartmentalModelMatrix *I_star;
            InitData *A0;
            double* p_rs;
            double *p_ir;
    };

    class FC_Beta : public FullConditional 
    {
        public:
            FC_Beta(CompartmentalModelMatrix *E_star, 
                    CompartmentalModelMatrix *S_star, 
                    InitData *A0,
                    CovariateMatrix *X,
                    double *p_se, 
                    double *beta,
                    double *rho); 

            CompartmentalModelMatrix *E_star; 
            CompartmentalModelMatrix *S_star; 
            InitData *A0;
            CovariateMatrix *X;
            double *p_se;
            double *beta;
            double *rho;
    };

    class FC_Rho : public FullConditional
    {
        public:
            FC_Rho(CompartmentalModelMatrix *S_star, 
                   CompartmentalModelMatrix *E_star, 
                   InitData *A0,
                   CovariateMatrix *X,
                   double *p_se, 
                   double *beta, 
                   double *rho
                   );
            CompartmentalModelMatrix *S_star; 
            CompartmentalModelMatrix *E_star; 
            InitData *A0;
            CovariateMatrix *X;
            double *p_se;
            double *beta 
            double *rho;
    }
    
    class FC_P_RS : public FullConditional
    {
        public:
            FC_P_RS(CompartmentalModelMatrix *S_star,
                    CompartmentalModelMatrix *R_star,
                    InitData *A0,
                    double *p_rs 
                    );
            CompartmentalModelMatrix *S_star;
            CompartmentalModelMatrix *R_star;
            InitData *A0;
            double *p_rs;

    };

    class FC_P_EI : public FullConditonal
    {
        public:
            FC_P_EI(CompartmentalModelMatrix *I_star,
                    CompartmentalModelMatrix *E_star,
                    InitData *A0,
                    double *p_ei);
            CompartmentalModelMatrix *I_star;
            CompartmentalModelMatrix *E_star;
            InitData *A0;
            double *p_ei;
    }

    class FC_P_IR : public FullConditional 
    {
        
        public:
            FC_P_IR(CompartmentalModelMatrix *I_star, 
                    CompartmentalModelMatrix *R_star,
                    InitData *A0,
                    double *p_ir);
            CompartmentalModelMatrix *I_star;
            CompartmentalModelMatrix *R_star;
            InitData *A0;
            double *p_ir;
    }

}
