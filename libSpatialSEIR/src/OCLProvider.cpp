#define __CL_ENABLE_EXCEPTIONS

#ifndef SPATIALSEIR_INCLUDEFILES
#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#endif


#include <CL/cl.hpp>
#include <OCLProvider.hpp>

// Implement OCLProvider Class
// Mostly placeholder code, based on Scarpino (2012)
SpatialSEIR::OCLProvider::OCLProvider()
{
    platforms = new std::vector<cl::Platform>();

    cl_uint i;
    try
    {
        cl::Platform::get(platforms);
        (*platforms)[0].getDevices(CL_DEVICE_TYPE_ALL, platformDevices);
        cl::Context context(*platformDevices);
        *ctxDevices = context.getInfo<CL_CONTEXT_DEVICES>();
        for (i = 0; i<ctxDevices -> size();i++)
        {
            deviceNames -> push_back((*ctxDevices)[i].getInfo<CL_DEVICE_NAME>());
            cout << "Adding Device: " 
                 << (*deviceNames)[i]
                 << endl;
        }
    }
    catch(cl::Error e)
    {
        cout << e.what() << ": Error Code "
             << e.err() << endl;
    }
}
