#define __CL_ENABLE_EXCEPTIONS

#ifndef SPATIALSEIR_INCLUDEFILES
#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#endif

#include <OCLProvider.hpp>

// Implement OCLProvider Class
// Mostly placeholder code, based on Scarpino (2012)
void SpatialSEIR::OCLProvider()
{
    platforms = new std::vector<cl::Platform>();

    cl_uint i;
    try
    {
        cl::Platform::get(&platforms);
        platforms[0].getDevices(CL_DEVICE_TYPE_ALL, platformDevices);
        cl::Context context(platformDevices);
        ctxDevices = context.getInfo<CL_CONTEXT_DEVICES>();
        for (i = 0; i<ctxDevices.size();i++)
        {
            deviceNames.push_back(ctxDevices[i]);
            cout << "Adding Device: " 
                 << deviceNames[i].c_str()
                 << endl;
        }
    }
    catch(cl::Errir e)
    {
        cout << e.what() << ": Error Code "
             << e.err() << endl;
    }
    return 0;
}
