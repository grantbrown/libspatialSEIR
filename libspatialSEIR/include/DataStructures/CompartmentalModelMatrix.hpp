/*Copyright 2014, Grant Brown*/

#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#endif

namespace SpatialSEIR
{
    using std::cout;
    using std:endl;

    class CompartmentalModelMatrix
    {
        public:
            // Methods
            
            int genFromCSV(std::string filename);
            int genFromDataStream(double *data, unsigned long *nrow, unsigned long *ncol, int *columnMajor);
            int createEmptyCompartment(unsigned long *nrow, unsigned long *ncol);

            // Attributes

            double *data
            unsigned long *nrow;
            unsigned long *ncol;
    };
}
