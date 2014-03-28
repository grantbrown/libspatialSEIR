/*Copyright 2014, Grant Brown*/

#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#endif

#ifndef COMPARTMENTAL_MODEL_MATRIX_INC
#define COMPARTMENTAL_MODEL_MATRIX_INC

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class CompartmentalModelMatrix
    {
        public:
            // Methods
            
            int genFromText(std::string filename);
            int genFromDataStream(int *indata, int *inrow, int *incol);
            int createEmptyCompartment(int *inrow, int *incol);
            ~CompartmentalModelMatrix();

            // Attributes

            int *data;
            int *nrow;
            int *ncol;

        private:
            int readDataFile(const char fn[]);
    };
}

#endif
