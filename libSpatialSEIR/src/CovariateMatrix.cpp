/*Copyright 2014, Grant Brown*/

#include <CovariateMatrix.hpp>
#include <cmath>
#ifdef LSS_USE_BLAS
	#include <cblas.h>
#endif
#include <LSS_FullConditional.hpp>
#include <IOProvider.hpp>



namespace SpatialSEIR
{
    int CovariateMatrix::genFromDataStream(double *indata_x,
                                           int inrow_x, 
                                           int incol_x)
    {
        int numToAlloc_x = (incol_x)*(inrow_x);

        X = new double[numToAlloc_x];
        int numVariables = (incol_x);
        decorrelationProjectionMatrix = new double[numVariables*numVariables];

        nrow_x = new int;
        ncol_x = new int;

        (*nrow_x) = (inrow_x);
        (*ncol_x) = (incol_x);
        int i; 
        for (i = 0; i < numToAlloc_x; i++)
        {
            X[i] = indata_x[i];
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
        int numVariables = (*ncol_x);
        int matrixRows;

        matrixRows = *nrow_x; 
        CovariateMatrix::MatrixMapType Amap(X, matrixRows, numVariables);
        CovariateMatrix::MatrixMapType outMap(decorrelationProjectionMatrix, numVariables, numVariables);

        CovariateMatrix::MatrixType Bmat = ((Amap.transpose() * Amap));
        CovariateMatrix::MatrixType Cmat = (Bmat.transpose() * Bmat);
        CovariateMatrix::MatrixType Dmat = (Bmat * Cmat.inverse() * Bmat.transpose()); 
        CovariateMatrix::MatrixType Emat = Eigen::MatrixXd::Identity(Dmat.rows(), Dmat.cols()) - Dmat;         
        outMap.noalias() = Emat;
    }

    int CovariateMatrix::calculate_eta_CPU(double *eta, double *beta)
    {
        try
        {
#ifdef LSS_USE_BLAS
            cblas_dgemv(CblasColMajor,
                        CblasNoTrans,
                        *nrow_x,
                        *ncol_x,
                        1.0,
                        X,
                        *nrow_x,
                        beta,
                        1,
                        0.0,
                        eta,
                        1);
#else
            CovariateMatrix::MatrixMapType Xmap(X, *(nrow_x), *(ncol_x));
            CovariateMatrix::MatrixMapType betaMap(beta, *ncol_x, 1);
            CovariateMatrix::MatrixMapType etaMap(eta, (*nrow_x), 1);
            etaMap.noalias() = Xmap*betaMap;
#endif
        }
        catch(int e)
        {
            lssCout << "Error calculating eta, code: " << e << "\n";
            return(-1);
        }
        return(0);
    }

    int CovariateMatrix::calculate_eta_OCL(double *eta, double *beta)
    {
        // Not implemented
        return(calculate_eta_CPU(eta, beta)); 
    }

    CovariateMatrix::~CovariateMatrix()
    {
        delete[] X;
        delete[] decorrelationProjectionMatrix;
        delete nrow_x;
        delete ncol_x;
    }

}
