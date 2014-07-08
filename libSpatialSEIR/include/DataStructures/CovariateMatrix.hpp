/*Copyright 2014, Grant Brown*/

#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#endif

#ifndef COVARIATE_MATRIX_INC
#define COVARIATE_MATRIX_INC

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    struct covariateArgs
    {
        double* inData_x;
        double* inData_z;
        int* inRow_x;
        int* inCol_x;
        int* inRow_z;
        int* inCol_z;
    };

    class CovariateMatrix
    {
        public:
            // Methods
        
            int genFromDataStream(double *indata_x, double *indata_z, 
                                  int *inrow_x, int *incol_x,
                                  int *inrow_z, int *incol_z);

            // Eta function for X only covariate structures. 
            int calculate_fixed_eta_CPU(double *eta, double *beta);
            // Eta functions for combined (beta, gamma) vectors.
            int calculate_eta_CPU(double *eta, double *beta);
            int calculate_eta_OCL(double *eta, double *beta); 

            // Eta functions for separate (beta), (gamma) vectors. 
            int calculate_eta_CPU(double *eta, double *beta, double *gamma);
            int calculate_eta_OCL(double *eta, double *beta, double *gamma); 
            ~CovariateMatrix();
            // Attributes
            double *X; // Time invariant covariates
            double *Z; // Time varying covariates
            double *offset;
            int *offsetLength;
            std::vector<std::string>* varnames;
            int *nrow_x;
            int *ncol_x;
            int *nrow_z;
            int *ncol_z;
            
    };
}

#endif
