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

    class ModelContext
    {
        public:
            //Methods
            ModelContext();            
            ModelContext(CompartmentalModelMatrix *I_star, 
                         CovariateMatrix *X);
            ~ModelContext(); 
        
            //Attributes
            OCLProvider *oclProvider; 
            FC_S_Star *S_star_fc;
            FC_E_Star *E_star_fc;
            FC_R_Star *R_star_fc;
            FC_Beta *beta_fc;
            FC_P_Rho *rho_fc;
            FC_P_RS *p_rs_fc;
            FC_P_EI *p_ei_fc;
            FC_P_IR *p_ir_fc;

    };
}
