#define __CL_ENABLE_EXCEPTIONS

#include <CL/cl.hpp>
#include <OCLProvider.hpp>
#include<iostream>
#include<fstream>
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
        programs = new std::vector<cl::Program>();
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
        context = new cl::Context(*platformDevices);
        *ctxDevices = context -> getInfo<CL_CONTEXT_DEVICES>();
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
        std::cout << "Attempting to compile test kernel:\n";
        buildProgramForKernel("test_kernel.cl", *ctxDevices);
    }
    catch(cl::Error e)
    {
        cout << e.what() << ": Error Code "
             << e.err() << endl;
    }
}


int SpatialSEIR::OCLProvider::buildProgramForKernel(std::string kernelFile, std::vector<cl::Device> devices)
{
    const char* progName = ( std::string(LSS_KERNEL_DIRECTORY).append(kernelFile)).c_str();
    std::cout << "Looking for kernel file here: " << progName << "\n";
    std::ifstream programFile(progName);
    std::string programString(std::istreambuf_iterator<char>(programFile), 
                             (std::istreambuf_iterator<char>()));
    std::cout << "Kernel: \n";
    std::cout << programString.c_str() << "\n";
    cl::Program::Sources source(1, std::make_pair(programString.c_str(), 
                                programString.length() + 1));
    cl::Program program(*context, source);
    program.build(devices);
    std::string log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
    std::cout << log << "\n";
    std::vector<cl::Kernel> kernels;
    program.createKernels(&kernels);
    unsigned int i;
    for (i = 0; i < kernels.size(); i++)
    {
        std::cout << kernels[i].getInfo<CL_KERNEL_FUNCTION_NAME>() << "\n";
    }
    return(0);
}


SpatialSEIR::OCLProvider::~OCLProvider()
{
    delete[] platforms;
    delete[] platformDevices;
    delete[] allDevices;
    delete[] ctxDevices;
    delete[] deviceNames;
    delete[] workSizes;
    delete[] doublePrecision;
    delete[] programs;
    delete context;
}


