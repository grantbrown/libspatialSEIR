/*Copyright 2014, Grant Brown*/

#include <cmath>
#include <DistanceMatrix.hpp>
#include <IOProvider.hpp>

namespace SpatialSEIR
{
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
        data = new double[numToAlloc];
        numLocations = new int;
        (*numLocations) = (*indim);
        int i; 
        for (i = 0; i < numToAlloc; i++)
        {
            data[i] = indata[i];
        }
        return(0);
    }
    int DistanceMatrix::scaledInvFunc_CPU(double phi)
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
                data[i] = std::sqrt(((tmp < phi && std::abs(tmp) > 0.1) ? 1/tmp : 0));
            }
        }
        return(0);
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
                data[i] = std::sqrt((tmp < phi && std::abs(tmp) > 0.1) ? 1/tmp : 0);
            }
        }
        return(0);
    }

    void DistanceMatrix::makeRowStochastic()
    {
        if (*hasAlloc!=1)
        {
            lssCout << "Distance Matrix: no data to standardize.\n";
            throw(-1);
        }
        int i,j;
        double rowTotal;
        for (i = 0; i < *numLocations; i++)
        {
            rowTotal = 0.0;
            for(j = 0; j< *numLocations; j++)
            {
                rowTotal += data[i + j*(*numLocations)];
            }
            for(j = 0; j< *numLocations; j++)
            {
                data[i + j*(*numLocations)] /=  rowTotal;
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
