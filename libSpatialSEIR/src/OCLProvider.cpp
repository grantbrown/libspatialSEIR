#define __CL_ENABLE_EXCEPTIONS

#include <CL/cl.hpp>
#include <OCLProvider.hpp>

#ifndef SPATIALSEIR_INCLUDEFILES
#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#endif


// Implement OCLProvider Class
// Mostly placeholder code, based on Scarpino (2012)
SpatialSEIR::OCLProvider::OCLProvider()
{
    try
    {
        platforms = new std::vector<cl::Platform>();
        platformDevices = new std::vector<cl::Device>();
        allDevices = new std::vector<cl::Device>();
        ctxDevices = new std::vector<cl::Device>();
        deviceNames = new std::vector<std::string>();
    }
    catch(cl::Error e)
    {
        cout << "Problem getting platforms:" << endl;
        cout << e.what() << ": Error Code " << e.err() << endl;
        return;
    }

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



SpatialSEIR::OCLProvider::~OCLProvider()
{
    delete[] platforms;
    delete[] platformDevices;
    delete[] allDevices;
    delete[] ctxDevices;
    delete[] deviceNames;
}


