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

    class CovariateMatrix
    {
        public:
            // Methods
            
            int genFromText(std::string filename);
            int genFromDataStream(double *indata, int *inrow, 
                    int *incol);
            ~CovariateMatrix();

            // Attributes
            double *data;
            std::vector<std::string>* varnames;
            int *nrow;
            int *ncol;

        private:
            int readDataFile(const char fn[]);

            

    };
}
