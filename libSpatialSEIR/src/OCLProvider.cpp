#define __CL_ENABLE_EXCEPTIONS

#include <CL/cl.hpp>
#include <OCLProvider.hpp>
#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>


// Implement OCLProvider Class
// Mostly placeholder code, based on Scarpino (2012)
SpatialSEIR::OCLProvider::OCLProvider()
{
    std::cout << "Setting up OpenCL Interface\n";
    try
    {
        platforms = new std::vector<cl::Platform>();
        platformDevices = new std::vector<cl::Device>();
        allDevices = new std::vector<cl::Device>();
        ctxDevices = new std::vector<cl::Device>();
        deviceNames = new std::vector<std::string>();
        workSizes = new std::vector<std::vector<size_t> >();
        doublePrecision = new std::vector<cl_uint>();
    }
    catch(cl::Error e)
    {
        cout << "Problem getting platforms:" << endl;
        cout << e.what() << ": Error Code " << e.err() << endl;
        return;
    }

    cl_uint i;
    cl_uint j;
    std::string clExt;
    try
    {
        cl::Platform::get(platforms);
        (*platforms)[0].getDevices(CL_DEVICE_TYPE_ALL, platformDevices);
        cl::Context context(*platformDevices);
        *ctxDevices = context.getInfo<CL_CONTEXT_DEVICES>();
        for (i = 0; i<ctxDevices -> size();i++)
        {
            deviceNames -> push_back((*ctxDevices)[i].getInfo<CL_DEVICE_NAME>());
            workSizes -> push_back((*ctxDevices)[i].getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>());
            clExt = (*ctxDevices)[i].getInfo<CL_DEVICE_EXTENSIONS>();
            doublePrecision -> push_back((clExt.find("cl_khr_fp64") != std::string::npos));
            cout << "   Adding Device: " 
                 << (*deviceNames)[i] << endl;
            cout << "      Workgroup Sizes: ";
            for (j = 0; j < (workSizes[i].size()); j++)
            {
                cout << ((*workSizes)[i])[j] << ", ";
            }
            cout << endl;
            cout << "      Supports Double Precision: " << (*doublePrecision)[i] << endl;
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


