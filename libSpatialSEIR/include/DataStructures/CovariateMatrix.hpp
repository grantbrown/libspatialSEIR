/*Copyright 2014, Grant Brown*/

#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#endif

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class CovariateMatrix
    {
        public:
            // Methods
        
            int genFromDataStream(double *indata_x, double *indata_z, 
                                  int *inrow_x, int *incol_x,
                                  int *inrow_z, int *incol_z);

            int calculate_eta_CPU(double *eta, double *beta, double *gamma);
            int calculate_eta_OCL(double *eta, double *beta, double *gamma); 
            ~CovariateMatrix();
            // Attributes
            double *X; // Time invariant covariates
            double *Z; // Time varying covariates
            std::vector<std::string>* varnames;
            int *nrow_x;
            int *ncol_x;
            int *nrow_z;
            int *ncol_z;
    };
}
