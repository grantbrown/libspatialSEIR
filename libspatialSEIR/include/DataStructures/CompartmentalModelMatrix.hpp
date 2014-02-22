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
    using std::endl;

    class CompartmentalModelMatrix
    {
        public:
            // Methods
            
            int genFromText(std::string filename);
            int genFromDataStream(double *indata, unsigned long *inrow, unsigned long *incol, int *columnMajor);
            int createEmptyCompartment(unsigned long *inrow, unsigned long *incol);
            ~CompartmentalModelMatrix();

            // Attributes

            int *data;
            unsigned long *nrow;
            unsigned long *ncol;

        private:
            int readDataFile(const char fn[]);
    };
}
