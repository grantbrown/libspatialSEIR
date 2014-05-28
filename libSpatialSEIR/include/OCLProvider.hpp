#ifndef OCL_PROVIDER_INC
#define OCL_PROVIDER_INC

#define __CL_ENABLE_EXCEPTIONS

#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#endif

#include <CL/cl.hpp>


namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class OCLProvider
    {
        public:  
            //Methods 
            OCLProvider();

            int buildProgramForKernel(std::string kernelFile, 
                    std::vector<cl::Device> devices);

            ~OCLProvider();
            //Attributes
            cl::Context *context;
            std::vector<cl::Platform> *platforms;
            std::vector<cl::Device> *platformDevices, 
                                    *allDevices, 
                                    *ctxDevices; 
            std::vector<std::string> *deviceNames;
            std::vector<std::vector<size_t> > *workSizes;
            std::vector<cl_uint> *doublePrecision;
            std::vector<cl::Program> *programs;
    };
}

#endif
