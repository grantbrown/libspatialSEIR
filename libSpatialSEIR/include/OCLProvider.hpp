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
            ~OCLProvider();
            //Attributes
            std::vector<cl::Platform> *platforms;
            std::vector<cl::Device> *platformDevices, 
                                    *allDevices, 
                                    *ctxDevices;
            std::vector<std::string> *deviceNames;
    };
}

#endif
