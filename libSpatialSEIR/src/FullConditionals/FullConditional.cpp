#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<FullConditional.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    /*
     * Helper functions
     */
    int matMult(double* output, double * A, double * B, int Arow, int Acol, 
            int Brow, int Bcol, bool TransA, bool TransB, int ldA, int ldB, int ldC)
    {
        // Use BLAS to matrix multiply, assume column major, non transposing.
        // Double check this code when I have internet access. 

        if ((TransA ? Arow : Acol) != (TransB ? Bcol : Brow))
        {
            std::cerr << "Invalid Matrix Dimensions: " << std::endl;
            std::cerr << "A: " <<  Arow << " x " << Acol << std::endl;
            std::cerr << "B: " <<  Brow << " x " << Bcol << std::endl;
            std::cerr << "Transpose: " << TransA << ", " << TransB << std::endl;
            throw(-1);
        }
        cblas_dgemm(CblasColMajor,
                    TransA ? CblasTrans : CblasNoTrans,
                    TransB ? CblasTrans : CblasNoTrans,
                    Arow, Bcol, Brow,
                    1.0, 
                    A, ldA,  
                    B, ldB,  
                    0.0, output, ldC);
        return 0; 
    }


    /*
     *
     * Implement the data container class InitData
     *
     */    
 
 
    InitData::InitData(int *_S0, 
                       int *_E0,
                       int *_I0,
                       int *_R0,
                       int *nLoc)
    {
        this -> populate(_S0, _E0, _I0, _R0, 
                nLoc);
    }

    InitData::InitData()
    {
        // Do nothing
    }

    void InitData::populate(int *_S0, 
                       int *_E0,
                       int *_I0,
                       int *_R0,
                       int *nLoc
                       )
    {
        S0 = new int[*nLoc];
        E0 = new int[*nLoc];
        I0 = new int[*nLoc];
        R0 = new int[*nLoc];
        numLocations = new int;
        *numLocations = *nLoc;
        int i;
        for (i = 0; i < *nLoc; i++)
        {
            S0[i] = _S0[i]; 
            E0[i] = _E0[i];
            I0[i] = _I0[i];
            R0[i] = _R0[i];
        }
    }

    InitData::~InitData()
    {
        delete S0;
        delete E0;
        delete I0;
        delete R0;
        delete numLocations;
    }

    double FullConditional::acceptanceRatio()
    {
        return((*accepted*1.0)/(*samples));
    }

    /*
     *
     * Auto-tuning logic
     *
     */

    void FullConditional::updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange)
    {
        if (*samples == 0)
        {
            std::cout << "No sampling data. \n";
            return;
        }
        if (proportionChange <= 0 || proportionChange >= 1)
        {
            std::cerr << "Invalid Proportion: " << proportionChange << "\n";
            throw(-1);
        }
        double currentRatio = (this -> acceptanceRatio());
        if ((currentRatio > desiredRatio) && (std::abs(currentRatio - desiredRatio) > targetWidth))
        {
            (*sliceWidth)*=(1+proportionChange); 
        }
        else if ((currentRatio < desiredRatio) && (std::abs(currentRatio - desiredRatio) > targetWidth))
        {
            (*sliceWidth)*=(1-proportionChange);           
        }

        *samples = 0;
        *accepted = 0;
    }
}


