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
                                                   double width,
                                                   double* metropolisComponent1,
                                                   double* metropolisComponent2)
    {        
        // Not implemented
    }

    void CompartmentFullConditional::proposeUpdate(int* initCompartment,
                                                   int* indexList,
                                                   int indexListLength,
                                                   double width,
                                                   double* metropolisComponent1,
                                                   double* metropolisComponent2)
    {
        // Not implemented
    }

    void ParameterFullConditional::proposeUpdate(double* variable,
                                                 int varLen,
                                                 double* width,
                                                 double* metropolisComponent1,
                                                 double* metropolisComponent2)
    {
        // Not Implemented
    }

    void ParameterFullConditional::proposeUpdate(double* variable,
                                                 int varLen,
                                                 CovariateMatrix* priorMatrix,
                                                 double* width,
                                                 double* metropolisComponent1,
                                                 double* metropolisComponent2)
    {
        // Not Implemented
    }

    
    void InitCompartmentFullConditional::proposeUpdate(int* initCompartment,
                                                       int initCompartmentLength,
                                                       double width,
                                                       double* metropolisComponent1,
                                                       double* metropolisComponent2)
    {
        // Not implemented
    }

    void InitCompartmentFullConditional::proposeUpdate(int* initCompartment,
                                                       int* indexList,
                                                       int indexListLength,
                                                       double width,
                                                       double* metropolisComponent1,
                                                       double* metropolisComponent2)
    {
        // Not implemented
    }


    void HybridFullConditional::proposeUpdate(double* variable,
                                              int varLen,
                                              int* destCompartment,
                                              double* varWidth,
                                              double compWidth,
                                              double* metropolisComponent1,
                                              double* metropolisComponent2)
    {
        // Not Implemented
    }

    void HybridFullConditional::proposeUpdate(double* variable,
                                              int varLen,
                                              int* destCompartment,
                                              double* varWidth,
                                              double compWidth,
                                              CovariateMatrix* priorMatrix,
                                              double* metropolisComponent1,
                                              double* metropolisComponent2)
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
                                              int indexListLength,
                                              double* metropolisComponent1,
                                              double* metropolisComponent2)
    {
        // Not Implemented
    }
}
