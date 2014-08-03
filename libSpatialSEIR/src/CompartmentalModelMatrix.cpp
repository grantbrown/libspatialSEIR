/*Copyright 2014, Grant Brown*/
#include <iostream>
#include <fstream>
#include <stdint.h>
#include <CompartmentalModelMatrix.hpp>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    int CompartmentalModelMatrix::genFromDataStream(int *indata,
                                                    int *inrow, 
                                                    int *incol)
    {
        int numToAlloc = (*incol)*(*inrow);
        data = new int[numToAlloc];
        nrow = new int;
        ncol = new int;
        (*nrow) = (*inrow);
        (*ncol) = (*incol);
        int i; 
        for (i = 0; i < (*incol)*(*inrow); i++)
        {
            data[i] = indata[i];
        }
        dataMatrix = new IntMatrixMapType(data, *nrow, *ncol);
        return(0);
    }

    int CompartmentalModelMatrix::createEmptyCompartment(int *inrow, int *incol)
    {
        nrow = new int;
        ncol = new int;
        *nrow = *inrow;
        *ncol = *incol;

        data = new int[(*nrow) * (*ncol)];
        int i;
        for (i = 0; i < ((*nrow)*(*ncol)); i ++ )
        {
            data[i] = 0;
        }
        dataMatrix = new IntMatrixMapType(data, *nrow, *ncol);
        return(0);
    }
    long unsigned int CompartmentalModelMatrix::marginSum(int margin, int slice)
    {
        int64_t output = 0; 
        int startIdx, i;
        // Rowsum
        if (margin == 1)
        {
            startIdx = slice;
            for (i = 0; i < (*ncol); i++)
            {
                output += data[startIdx + i*(*nrow)];    
            }
            return(output);
        }
        // Colsum
        else if (margin == 2)
        {
            startIdx = slice*(*nrow);
            for (i = 0; i < (*nrow); i++)
            {
                output += data[startIdx + i];
            }
            return(output);
        }
        // Anything else indicates a total sum
        int maxIdx = (*nrow)*(*ncol);
        for (i = 0; i < maxIdx; i++)
        {
            output += data[i];
        }
        return(output);
    }

    CompartmentalModelMatrix::~CompartmentalModelMatrix()
    {
        delete dataMatrix;
        delete[] data;
        delete nrow;
        delete ncol;
    }

}
