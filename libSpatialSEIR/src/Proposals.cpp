#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<LSS_FullConditional.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    
    
    void CompartmentFullConditional::proposeUpdate(int* initCompartment,
                                                   int initCompartmentLength,
                                                   double width)
    {
        // Not implemented
    }

    void CompartmentFullConditional::proposeUpdate(int* initCompartment,
                                                   int* indexList,
                                                   int indexListLength,
                                                   double width)
    {
        // Not implemented
    }

    void ParameterFullConditional::proposeUpdate(double* variable,
                                                 int varLen,
                                                 double* width)
    {
        // Not Implemented
    }

    void ParameterFullConditional::proposeUpdate(double* variable,
                                                 int varLen,
                                                 CovariateMatrix* priorMatrix,
                                                 double* width)
    {
        // Not Implemented
    }

    
    void InitCompartmentFullConditional::proposeUpdate(int* initCompartment,
                                                       int initCompartmentLength,
                                                       double width)
    {
        // Not implemented
    }

    void InitCompartmentFullConditional::proposeUpdate(int* initCompartment,
                                                       int* indexList,
                                                       int indexListLength,
                                                       double width)
    {
        // Not implemented
    }


    void HybridFullConditional::proposeUpdate(double* variable,
                                              int varLen,
                                              int* destCompartment,
                                              double* varWidth,
                                              double compWidth)
    {
        // Not Implemented
    }

    void HybridFullConditional::proposeUpdate(double* variable,
                                              int varLen,
                                              int* destCompartment,
                                              double* varWidth,
                                              double compWidth,
                                              CovariateMatrix* priorMatrix)
    {
        // Not Implemented
    }

    void HybridFullConditional::proposeUpdate(double* variable,
                                              int varLen,
                                              int* destCompartment,
                                              double* varWidth,
                                              double compWidth,
                                              CovariateMatrix* priorMatrix,
                                              int *indexList,
                                              int indexListLength)
    {
        // Not Implemented
    }
}
