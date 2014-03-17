#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<cstring>
#include<vector>
#endif

#ifndef OCL_PROVIDER_INC
#define OCL_PROVIDER_INC
#include "OCLProvider.hpp"
#endif

#ifndef COMPARTMENTAL_MODEL_MATRIX_INC
#define COMPARTMENTAL_MODEL_MATRIX_INC
#include "CompartmentalModelMatrix.hpp"
#endif

#ifndef COVARIATE_MATRIX_INC
#define COVARIATE_MATRIX_INC
#include "CovariateMatrix.hpp"
#endif

#ifndef FULL_CONDITIONAL_INC
#define FULL_CONDITIONAL_INC
#include "FullConditional.hpp"
#endif

#ifndef DISTANCEMATRIX
#define DISTANCEMATRIX 
#include "DistanceMatrix.hpp"
#endif

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class ModelContext
    {
        public:
            //Methods
            ModelContext();
            ~ModelContext(); 

            // Method: calculateS
            // Accesses: A0, S_star, E_star
            // Updates: S
            void calculateS_CPU();
            void calculateS_OCL();

            // Method: calculateE
            // Accesses: A0, I_star, E_star
            // Updates: E
            void calculateE_CPU();
            void calculateE_OCL();

            // Method: calculateI
            // Accesses: A0, I_star, R_star
            // Updates: I
            void calculateI_CPU();
            void calculateI_OCL();

            // Method: calculateR
            // Accesses: A0, R_star, S_star
            // Updates: R
            void calculateR_CPU();
            void calculateR_OCL();

            // Method: calculatePi
            // Accesses: beta, I, N, distMat, rho
            // Updates: p_se
            void calculateP_SE_CPU();
            void calculateP_SE_OCL();
        
            //Logic provider classes
            OCLProvider *oclProvider; 
            FC_S_Star *S_star_fc;
            FC_E_Star *E_star_fc;
            FC_R_Star *R_star_fc;
            FC_Beta *beta_fc;
            FC_Rho *rho_fc;
            FC_P_RS *p_rs_fc;
            FC_P_EI *p_ei_fc;
            FC_P_IR *p_ir_fc;

            //Data
            CompartmentalModelMatrix* S;
            CompartmentalModelMatrix* E;
            CompartmentalModelMatrix* I;
            CompartmentalModelMatrix* R;
            CompartmentalModelMatrix* S_star;
            CompartmentalModelMatrix* E_star;
            CompartmentalModelMatrix* I_star;
            CompartmentalModelMatrix* R_star;
            InitData* A0;
            CovariateMatrix* X;
            DistanceMatrix* rawDistMat;
            DistanceMatrix* scaledDistMat;
            double* p_se;
            double* p_ei;
            double* p_ir;
            double* p_rs;
            int* N;// Matrix?
    };
}
