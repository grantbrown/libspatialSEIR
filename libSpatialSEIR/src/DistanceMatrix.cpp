/*Copyright 2014, Grant Brown*/


#include <iostream>
#include <fstream>
#include <cmath>

#ifndef DISTANCEMATRIX 
#define DISTANCEMATRIX
#include <DistanceMatrix.hpp>
#endif

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    DistanceMatrix::DistanceMatrix()
    {
        this -> hasAlloc = new int;
        *hasAlloc = 0;
    }

    int DistanceMatrix::genFromDataStream(double *indata,
                                          int *indim)
    {
        *hasAlloc = 1;
        int numToAlloc = (*indim)*(*indim);
        this -> data = new double[numToAlloc];
        this -> numLocations = new int;
        (*numLocations) = (*indim);
        int i; 
        for (i = 0; i < numToAlloc; i++)
        {
            data[i] = indata[i];
        }
    }
    int DistanceMatrix::scaledInvFunc_CPU(double phi, double* indata)
    {
        if (*hasAlloc != 1)
        {
            // Data not allocated, not enough info to allocate
            throw(-1);
        }
        else
        {
            int i; int maxVal = (*numLocations)*(*numLocations);
            double tmp;
            for (i = 0; i < maxVal; i++)
            {
                tmp = data[i];
                data[i] = ((tmp < phi && std::abs(tmp) > 0.1) ? 1/tmp : 0);
            }
        }
    }
    int DistanceMatrix::scaledInvFunc_CPU(double phi, double* indata, int* indim)
    {
        if (*hasAlloc != 1)
        {
            this -> genFromDataStream(&*indata, &*indim);
        }
        else
        {
            int i; int maxVal = (*indim)*(*indim);
            double tmp;
            for (i = 0; i < maxVal; i++)
            {
                tmp = data[i];
                data[i] = ((tmp < phi && std::abs(tmp) > 0.1) ? 1/tmp : 0);
            }
        }
    }

    DistanceMatrix::~DistanceMatrix()
    {
        if ((*hasAlloc) == 1)
        {
            delete[] data;
            delete numLocations;
        }
        delete hasAlloc;
    }

}
