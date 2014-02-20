#define __CL_ENABLE_EXCEPTIONS

#include <iostream>
#include <CL/cl.hpp>
#include <vector>

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
