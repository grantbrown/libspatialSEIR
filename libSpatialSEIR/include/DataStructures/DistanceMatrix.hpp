/*Copyright 2014, Grant Brown*/

#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#endif

#ifndef DISTANCE_MATRIX_INC
#define DISTANCE_MATRIX_INC

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    struct distanceArgs
    {
        double* inData;
        int* dim; 
    };

    struct scaledDistanceArgs
    {
        double* phi;
        double* inData;
        int* dim;
    };

    // Begin with a naive implementation, this should eventually 
    // be updated to take advantage of the inherent sparsity of 
    // a distance matrix with a max range, as well as the structure 
    // given by symmetry.
    //
    // First priority: make it work
    class DistanceMatrix
    {
        public:
            // Methods 
            DistanceMatrix();

            int genFromDataStream(double *indata, int *dim); 

            // Scaled options for already allocated distance matrix
            int scaledInvFunc_CPU(double phi, double *indata);
            int scaledInvFunc_OCL(double phi, double *indata);
            int gravityFunc_CPU(double phi, double *indata);
            int gravityFunc_OCL(double phi, double *indata);

            // Scaled options for unallocated distance matrix
            int scaledInvFunc_CPU(double phi, double *indata, int *dim);
            int scaledInvFunc_OCL(double phi, double *indata, int *dim);
            int gravityFunc_CPU(double phi, double *indata, int *dim);
            int gravityFunc_OCL(double phi, double *indata, int *dim);


            ~DistanceMatrix();
            // Attributes
            double *data;
            int *numLocations;

        private:
            int *hasAlloc;;
    };
}

#endif
