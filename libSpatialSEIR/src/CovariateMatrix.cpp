/*Copyright 2014, Grant Brown*/

#include <CovariateMatrix.hpp>
#include <iostream>
#include <sstream>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    int CovariateMatrix::genFromDataStream(double *indata_x,
                                           double *indata_z,
                                           int *inrow_x, 
                                           int *incol_x,
                                           int *inrow_z,
                                           int *incol_z)
    {
        int numToAlloc_x = (*incol_x)*(*inrow_x);
        int numToAlloc_z = (*incol_z)*(*inrow_z);

        double* X = new double[numToAlloc_x];
        double* Z = new double[numToAlloc_z];

        double* nrow_x = new double;
        double* ncol_x = new double;
        double* nrow_z = new double;
        double* ncol_z = new double;

        (*nrow_x) = (*inrow_x);
        (*ncol_x) = (*incol_x);
        (*nrow_z) = (*inrow_z);
        (*ncol_z) = (*incol_z);

        int i; 
        for (i = 0; i < numToAlloc_x; i++)
        {
            X[i] = indata_x[i];
        }
        for (i = 0; i < numToAlloc_z; i++)
        {
            Z[i] = indata_z[i];
        }

    }

    CovariateMatrix::~CovariateMatrix()
    {
        std::cout << "Covariate Matrix Destroyed" << std::endl;
        delete[] X;
        delete[] Z;
        delete[] nrow_x;
        delete[] ncol_x;
        delete[] nrow_z;
        delete[] ncol_z;

    }

}
