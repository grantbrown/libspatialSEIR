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
        cpuQueue = new cl::CommandQueue();
        gpuQueue = new cl::CommandQueue();
        doublePrecision = new std::vector<cl_uint>();
        programs = new std::vector<cl::Program>();
        test_kernel = new cl::Kernel();
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
        // This might be a memory leak of the vector containing test kernel, though low 
        // priority as it runs once and the bulk of the data is tracked (via the pointer 
        // test_kernel)
        //
        // todo: figure out more about whether or not this stuff would work as pass by value, 
        // and how much OpenCL boilerplate needs to remain in memory for a particular kernel 
        // to be useful (ie, can we discard program objects after obtaining kernels?). Most 
        // tutorials and documentation don't go into this issue, so we probably need to dig into 
        // cl.hpp to verify. 
        test_kernel = &((*(buildProgramForKernel("test_kernel.cl", *ctxDevices)))[0]);

        *cpuQueue = cl::CommandQueue(*context, (*ctxDevices)[0]); 

        test();
    }
    catch(cl::Error e)
    {
        cout << e.what() << ": Error Code "
             << e.err() << endl;
    }
}


std::vector<cl::Kernel>* SpatialSEIR::OCLProvider::buildProgramForKernel(std::string kernelFile, std::vector<cl::Device> devices)
{
    int err;
    const char* progName = ( std::string(LSS_KERNEL_DIRECTORY).append(kernelFile)).c_str();
    std::ifstream programFile(progName);
    std::string programString(std::istreambuf_iterator<char>(programFile), 
                             (std::istreambuf_iterator<char>()));
    cl::Program::Sources source(1, std::make_pair(programString.c_str(), 
                                programString.length() + 1));
    cl::Program* program = new cl::Program(*context, source);
    err = program -> build(devices);
    std::string log = program -> getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
    if (err != 0)
    {
        std::cerr << "Error building OpenCL Kernel, code: " << err << "\n"; 
        std::cerr << "Build Log: \n" << log << "\n";
        std::cerr << "Kernel Source: \n" << programString.c_str() << "\n";
        throw(-1);
    }
    std::vector<cl::Kernel>* kernels = new std::vector<cl::Kernel>();
    program -> createKernels(kernels);
    programs -> push_back(*program);
    return(kernels);
}

int SpatialSEIR::OCLProvider::test()
{

    float* A = new float[100];
    float* B = new float[100];
    float* out = new float[100];
    float* cpuOut = new float[100];
    double bias=0.0;
    int err;
    void* mappedMemory;
    int i;
    for (i = 0; i < 100; i++)
    {
        A[i] = i*1.0;
        B[i] = i*1.0;
        cpuOut[i] = A[i] + B[i];
    }

    cl::Buffer bufferA(*context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 100*sizeof(float), A);
    cl::Buffer bufferB(*context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 100*sizeof(float), B);
    cl::Buffer bufferOut(*context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 100*sizeof(float), out);

    err = test_kernel -> setArg(0, bufferA);
    err |= test_kernel -> setArg(1, bufferB);
    err |= test_kernel -> setArg(2, bufferOut);
    if (err < 0)
    {
        std::cerr << "Couldn't set kernel args.\n"; 
        throw(err);
    }

    cpuQueue -> enqueueNDRangeKernel(*test_kernel, 
                                   0,           // Global Offset
                                   100,         // Global Size
                                   1,        // Local Size
                                   NULL,        // Event Vector
                                   NULL         // Event Pointer
                                   );

    mappedMemory = (cpuQueue -> enqueueMapBuffer(bufferOut, CL_TRUE, CL_MAP_READ, 0, 100*sizeof(float)));

    memcpy(out, mappedMemory, sizeof(float)*100);

    for (i=0;i<100;i++)
    {
        bias += pow(cpuOut[i] - out[i], 2);
    }
    if (bias > 1e-20)
    {
        std::cerr << "OPEN CL ERROR: Invalid Test Result, bias = " << bias << "\n";
        throw(-1);
    } 
    // todo: test output
    cpuQueue -> enqueueUnmapMemObject(bufferOut, mappedMemory);

    delete A;
    delete B;
    delete out;
    delete cpuOut;
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
    delete test_kernel;
    delete context;
    delete cpuQueue;
    delete gpuQueue;
}


