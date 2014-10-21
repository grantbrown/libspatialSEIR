/*Copyright 2014, Grant Brown*/

#include <CovariateMatrix.hpp>
#include <cmath>
#ifdef LSS_USE_BLAS
	#include <cblas.h>
#endif
#include <LSS_FullConditional.hpp>
#include <Eigen/Core>
#include <Eigen/LU>
#include <IOProvider.hpp>

namespace SpatialSEIR
{
    int CovariateMatrix::genFromDataStream(double *indata_x,
                                           double *indata_z, 
                                           int *inrow_x, 
                                           int *incol_x,
                                           int *inrow_z,
                                           int *incol_z)
    {
        int numToAlloc_x = (*incol_x)*(*inrow_x);
        int numToAlloc_z;
        if (*incol_z == 0)
        {
            numToAlloc_z = (*inrow_z);
        }
        else
        {
            numToAlloc_z = (*incol_z)*(*inrow_z);
        }


        X = new double[numToAlloc_x];
        Z = new double[numToAlloc_z];
        int numVariables = (*incol_x)+(*incol_z);
        decorrelationProjectionMatrix = new double[numVariables*numVariables];

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
        buildDecorrelationProjectionMatrix();
        return(0);
    }

    void CovariateMatrix::buildDecorrelationProjectionMatrix()
    {
        // There's no use parallelizing this. Even though it might be 
        // computationally intensive, this operation only needs to be done once per 
        // chain per covariate matrix. 

        // Combine X and Z appropriately to one large covariate matrix X
        // We want a projection onto N(X'X)
        //
        // The projection onto the column space of A= X'X = P(A) = (A) %*% (A'A)^(-1)) %*% (A')

        // The projection onto the null space of X' is then I - P(A)

        // The steps to compute I-P(A) are therefore:
        // 1. Calculate B = (A'A)
        // 2. Calculate C = B^(-1)
        // 3. Calculate D= A'CA
        // 5. Calculate E = I - D


        // Step 0. Create big matrix X.  
        int numVariables = (*ncol_x)+(*ncol_z);
        double *bigX;
        int matrixRows;
        int i;
        // Case 1: X matrix only, time varying. 
        if (*nrow_z == 0)
        {
            matrixRows = *nrow_x;
            bigX = new double[numVariables*(matrixRows)];
            memcpy(bigX, X, numVariables*(*nrow_x)*sizeof(double));    
        }
        // Case 2: X and Z matrices, fixed and time varying
        else
        {
            matrixRows = *nrow_z;
            bigX = new double[numVariables*(matrixRows)];
            memset(bigX, 0, numVariables*(*nrow_z)*sizeof(double));    

            int numLoc = *nrow_x;
            int numTpt = (*nrow_z)/numLoc;
            int xCol, xRow;
            int idx1, idx2;
            // Load up fixed covariates
            for (xCol = 0; xCol < (*ncol_x); xCol++)
            {
                for (xRow = 0; xRow < numLoc; xRow ++)
                {
                    for (i = 0; i < numTpt; i++)
                    {
                        idx1 = xCol*(*nrow_z) + xRow*numTpt + i;
                        idx2 = xCol*numLoc + xRow;
                        bigX[idx1] = X[idx2];
                    }
                }
            }
            // Load up time varying covariates
            int zCol, zRow;
            for (zCol = 0; zCol < (*ncol_z); zCol ++)
            {
                for (zRow = 0; zRow < (*nrow_z); zRow ++) 
                {
                    idx1 = (*ncol_x + zCol)*(*nrow_z) + zRow;
                    idx2 = zCol*(*nrow_z) + zRow;
                    bigX[idx1] = Z[idx2];
                }
            } 
        }

        // Step 1. Calculate A = (XX')
        double* A = bigX;
        
        typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatrixType;
        typedef Eigen::Map<MatrixType, Eigen::ColMajor> MatrixMapType;

        MatrixMapType Amap(A, matrixRows, numVariables);
        MatrixMapType outMap(decorrelationProjectionMatrix, numVariables, numVariables);


        MatrixType Bmat = ((Amap.transpose() * Amap));
        MatrixType Cmat = (Bmat.transpose() * Bmat);
        MatrixType Dmat = (Bmat * Cmat.inverse() * Bmat.transpose()); 
        MatrixType Emat = Eigen::MatrixXd::Identity(Dmat.rows(), Dmat.cols()) - Dmat;         
        outMap.noalias() = Emat;
 
        delete[] bigX;
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
                   eta[j] += X[j + i*(*nrow_x)]*beta[i];
               }
            }
        }
        catch(int e)
        {
            lssCout << "Error calculating eta, code: " << e << "\n";
            return(-1);
        }
        return(0);
    }

    // Combined beta/gamma parameters
    int CovariateMatrix::calculate_eta_CPU(double *eta, double *beta)
    {
        // TODO: do this in LAPACK
        // TODO: bite the bullet and just store this as one big matrix, 
        // forget the distinction between time varying and not. Storage is cheap, 
        // and the extra overhead is confusing. This will require API cleanup. 
        try
        {
            int nLoc = (*nrow_x);
            int nTpt = (*nrow_z)/(*nrow_x);
            int compIdx;
            int i; int j;
            int k;
            // Initialize eta 
            for (i = 0; i < (*nrow_z); i++)
            {
                eta[i] = 0.0;
            }
            // Locations
            for (i = 0; i < nLoc; i++)
            {
                compIdx = i*nTpt;
                // Time points
                for (j = 0; j < nTpt; j++)
                {
                    // Fixed Co-variates
                    for (k = 0; k < (*ncol_x); k++)
                    {
                        eta[compIdx] += X[i + k*nLoc]*beta[k];
                    }
                    // Time Varying Co-variates
                    for (k = 0; k < (*ncol_z); k++)
                    {
                        eta[compIdx] += Z[j + i*nTpt + k*(*nrow_z)]*beta[k + (*ncol_x)];
                    }
                    compIdx++;
                }
            }
        }
        catch(int e)
        {
            lssCout << "Error calculating eta, code: " << e << "\n";
            return(-1);
        }
        return(0);
    }

    CovariateMatrix::~CovariateMatrix()
    {
        delete[] X;
        delete[] Z;
        delete[] decorrelationProjectionMatrix;
        delete nrow_x;
        delete ncol_x;
        delete nrow_z;
        delete ncol_z;

    }

}
