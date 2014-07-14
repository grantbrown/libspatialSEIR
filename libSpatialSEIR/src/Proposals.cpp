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

    
    
    void CompartmentFullConditional::proposeUpdate(int* compartment,
                                                   int compartmentLength,
                                                   double width,
                                                   double* metropolisComponent1,
                                                   double* metropolisComponent2)
    {        
        int i;
        int x0, x1;
        for (i = 0; i < compartmentLength; i++)
        {
            x0 = compartment[i];
            x1 = std::floor(context -> random -> normal(x0 + 0.5, width));
            compartment[i] = x1;
            *metropolisComponent1 += (context -> random -> dnorm(x1, x0, width));
            *metropolisComponent2 += (context -> random -> dnorm(x0, x1, width));
        }
    }

    void CompartmentFullConditional::proposeUpdate(int* compartment,
                                                   int* indexList,
                                                   int indexListLength,
                                                   double width,
                                                   double* metropolisComponent1,
                                                   double* metropolisComponent2)
    {
        int i;
        int x0, x1;
        for (i = 0; i < indexListLength; i++)
        {
            x0 = compartment[indexList[i]];
            x1 = std::floor((*context) -> random -> normal(x0 + 0.5, width));
            compartment[indexList[i]] = x1;
            *metropolisComponent1 += (context -> random -> dnorm(x1, x0, width));
            *metropolisComponent2 += (context -> random -> dnorm(x0, x1, width));
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
}
