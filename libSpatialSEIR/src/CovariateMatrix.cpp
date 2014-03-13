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

        this -> X = new double[numToAlloc_x];
        this -> Z = new double[numToAlloc_z];

        this -> nrow_x = new int;
        this -> ncol_x = new int;
        this -> nrow_z = new int;
        this -> ncol_z = new int;

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

    int CovariateMatrix::calculate_eta_CPU(double *eta, double *beta, double *gamma)
    {
        // This is a naiive, but hopefully correct implementation. 
        // TODO: do this in LAPACK
        try
        {

            int i; int j;
            // Initialize eta 
            for (i = 0; i < ((*ncol_x) + (*ncol_z)); i++)
            {
                eta[i] = 0.0;
            }
            // Calc eta_x
            for (i = 0; i < (*ncol_x); i++)
            {
                for (j = 0; j < (*nrow_x); j++) 
                {
                    eta[i] += X[j + i*(*nrow_x)]*beta[i];                    
                }
            }

            // Calc eta_z
            for (i = 0; i < (*ncol_z); i++)
            {
                for (j = 0; j < (*nrow_z); j++)
                {
                    eta[i+(*ncol_x)] += Z[j + i*(*ncol_z)]*gamma[i]; 
                }
            }

        }
        catch(int e)
        {
            std::cout << "Error calculating eta, code: " << e << std::endl;
            return(-1);
        }
        return(0);
    }

    int CovariateMatrix::calculate_eta_OCL(double *eta, double *beta, double *gamma)
    {
        //Not implemented
        return(-1);
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
