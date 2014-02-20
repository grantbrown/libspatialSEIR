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

    class OCLProvider()
    {
        public:  
            //Methods 
            void OCLProvider();
            //Attributes
            vector<cl::Platform> platform;
            vector<cl::Device> platformDevices, allDevices, ctxDevices;
            vector<cl::string> device_names;
    }
}
