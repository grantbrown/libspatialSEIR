#define __CL_ENABLE_EXCEPTIONS

#include <CL/cl.hpp>
#include <OCLProvider.hpp>
#include<iostream>
#include<fstream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>


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
        globalMemSizes = new std::vector<size_t>();
        localMemSizes = new std::vector<cl_ulong>();
        numComputeUnits = new std::vector<size_t>();
        cpuQueue = new cl::CommandQueue();
        gpuQueue = new cl::CommandQueue();
        doublePrecision = new std::vector<cl_uint>();
        programs = new std::vector<cl::Program>();
        test_kernel = new cl::Kernel();
        R_Star_p1_kernel = new cl::Kernel();
        R_star_args = new FC_R_Star_KernelData();
        R_star_args -> totalWorkUnits = -1;
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
        // Todo: Handle multiple platforms. 
        (*platforms)[0].getDevices(CL_DEVICE_TYPE_ALL, platformDevices);
        context = new cl::Context(*platformDevices);
        *ctxDevices = context -> getInfo<CL_CONTEXT_DEVICES>();
        for (i = 0; i<ctxDevices -> size();i++)
        {
            deviceNames -> push_back((*ctxDevices)[i].getInfo<CL_DEVICE_NAME>());
            workSizes -> push_back((*ctxDevices)[i].getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>());
            globalMemSizes -> push_back((*ctxDevices)[i].getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>());
            localMemSizes -> push_back((*ctxDevices)[i].getInfo<CL_DEVICE_LOCAL_MEM_SIZE>());
            numComputeUnits -> push_back((*ctxDevices)[i].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>());
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
            cout << "      Global Memory: " << (*globalMemSizes)[globalMemSizes -> size() - 1] << endl;
            cout << "      Local Memory: " << (*localMemSizes)[localMemSizes -> size() - 1] << endl;
            cout << "      Compute Units: " << (*numComputeUnits)[numComputeUnits -> size() - 1] << endl;
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

        R_Star_p1_kernel = &((*(buildProgramForKernel("R_Star_FC_Part1.cl", *ctxDevices)))[0]);
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
    int err = 1;
    std::string log;
    std::string LKD(LSS_KERNEL_DIRECTORY);
    LKD = LKD.append(kernelFile);
    const char* progName = LKD.c_str();

    std::ifstream programFile(progName);
    std::string programString(std::istreambuf_iterator<char>(programFile), 
                             (std::istreambuf_iterator<char>()));
    cl::Program::Sources source(1, std::make_pair(programString.c_str(), 
                                programString.length() + 1));

    cl::Program* program = new cl::Program(*context, source);
    std::vector<cl::Kernel>* kernels = new std::vector<cl::Kernel>();
    try
    {
        err = program -> build(devices);
        log = program -> getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
        program -> createKernels(kernels);
        programs -> push_back(*program);
    }
    catch(cl::Error e)
    {
        std::cout << "CL Error in: " << e.what()<< "\n";
        std::cout << "CL Error: " << e.err()<< "\n";
        log = program -> getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
        err = e.err();
    }
    if (err != 0)
    {
        std::cerr << "Error building OpenCL Kernel, code: " << err << "\n"; 
        std::cerr << "Looking for kernel file here: " << progName << "\n";
        std::cerr << "Build Log: \n" << log << "\n";
        std::cerr << "Kernel Source: \n" << programString.c_str() << "\n";
        throw(-1);
    }

    return(kernels);
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
    delete[] globalMemSizes;
    delete[] localMemSizes;
    delete[] numComputeUnits;
    delete R_star_args;
    delete test_kernel;
    delete R_Star_p1_kernel;
    delete context;
    delete cpuQueue;
    delete gpuQueue;
}


