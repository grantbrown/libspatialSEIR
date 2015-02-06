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

#include <Eigen/Core>
#include <Eigen/LU>
namespace SpatialSEIR
{

    struct covariateArgs
    {
        double* inData_x;
        int inRow_x;
        int inCol_x;
    };

    class CovariateMatrix
    {
        typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatrixType;
        typedef Eigen::Map<MatrixType, Eigen::ColMajor> MatrixMapType;
        public:
            // Methods
        
            int genFromDataStream(double *indata_x, 
                                  int inrow_x, int incol_x);

            int calculate_eta_CPU(double *eta, double *beta);
            int calculate_eta_OCL(double *eta, double *beta); 
            void buildDecorrelationProjectionMatrix();

            ~CovariateMatrix();
            // Attributes
            double *X; // Time invariant covariates
            double *decorrelationProjectionMatrix;
            double *offset;
            int *offsetLength;
            std::vector<std::string>* varnames;
            int *nrow_x;
            int *ncol_x; 
    };
}

#endif
