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

        X = new double[numToAlloc_x];
        Z = new double[numToAlloc_z];

        nrow_x = new int;
        ncol_x = new int;
        nrow_z = new int;
        ncol_z = new int;


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
        return(0);
    }

    int CovariateMatrix::calculate_fixed_eta_CPU(double *eta, double *beta)
    {
        // This is a naiive, but hopefully correct implementation. 
        // TODO: do this in LAPACK
        try
        {
            int i; int j;
            // Initialize eta 
            for (i = 0; i < (*nrow_x); i++)
            {
                eta[i] = 0.0;
            }
            for (j = 0; j < (*nrow_x); j++)
            {
               for (i = 0; i < (*ncol_x); i++) 
               {
                   eta[j] += X[j%(*nrow_x) + i*(*nrow_x)]*beta[i];
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

    // Combined beta/gamma parameters
    int CovariateMatrix::calculate_eta_CPU(double *eta, double *beta)
    {
        // This is a naiive, but hopefully correct implementation. 
        // TODO: do this in LAPACK
        try
        {
            int i; int j;
            // Initialize eta 
            for (i = 0; i < (*nrow_z); i++)
            {
                eta[i] = 0.0;
            }
            for (j = 0; j < (*nrow_z); j++)
            {
               for (i = 0; i < (*ncol_x); i++) 
               {
                   eta[j] += X[j%(*nrow_x) + i*(*nrow_x)]*beta[i];
               }
               for (i = 0; i < (*ncol_z); i++)
               {
                   eta[j] += Z[j + i*(*nrow_z)]*beta[i + (*ncol_x)];
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

    int CovariateMatrix::calculate_eta_CPU(double *eta, double *beta, double *gamma)
    {
        // This is a naiive, but hopefully correct implementation. 
        // TODO: do this in LAPACK
        try
        {
            int i; int j;

            // Initialize eta 
            for (i = 0; i < (*nrow_z); i++)
            {
                eta[i] = 0.0;
            }
            for (j = 0; j < (*nrow_z); j++)
            {
               for (i = 0; i < (*ncol_x); i++) 
               {
                   eta[j] += X[j%(*nrow_x) + i*(*ncol_x)]*beta[i];
               }
               for (i = 0; i < (*ncol_z); i++)
               {
                   eta[j] += Z[j + i*(*nrow_z)]*gamma[i];
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
