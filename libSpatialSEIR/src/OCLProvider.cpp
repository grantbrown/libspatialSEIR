#define __CL_ENABLE_EXCEPTIONS

#include <CL/cl.hpp>
#include <OCLProvider.hpp>
#include<iostream>
#include<fstream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>


SpatialSEIR::DeviceContainer::DeviceContainer(cl::Device inDevice, cl::Context inContext)
{
    device = new cl::Device;
    *device = inDevice; 
    commandQueue = new cl::CommandQueue();
    *commandQueue = cl::CommandQueue(inContext, inDevice); 
}
SpatialSEIR::DeviceContainer::~DeviceContainer()
{
    delete commandQueue;
    delete device;
}

SpatialSEIR::PlatformContainer::PlatformContainer(cl::Platform inPlatform)
{
    //Store platform
    platform = new cl::Platform(); 
    *platform = inPlatform;

    //Temporarily store all platform devices to create context
    std::vector<cl::Device> allDevices;
    platform -> getDevices(CL_DEVICE_TYPE_CPU, &allDevices);
    context = new cl::Context(allDevices);

    // Process and store devices by type
    devices = new std::vector<DeviceContainer>();
    deviceTypes = new std::vector<std::string>();
    deviceNames = new std::vector<std::string>();
    doublePrecision = new std::vector<cl_uint>();

    std::vector<cl::Device> cpuDevices; 
    std::vector<cl::Device> gpuDevices; 

    platform -> getDevices(CL_DEVICE_TYPE_CPU, &cpuDevices);
    platform -> getDevices(CL_DEVICE_TYPE_GPU, &gpuDevices);

    std::string clExt;
    cl_uint i;
    for (i = 0; i < cpuDevices.size(); i++)
    {   
        devices -> push_back(DeviceContainer(cpuDevices[i], *context));
        clExt = cpuDevices[i].getInfo<CL_DEVICE_EXTENSIONS>();
        doublePrecision -> push_back((clExt.find("cl_khr_fp64") != std::string::npos));
        deviceTypes -> push_back("CPU");
        deviceNames -> push_back(cpuDevices[i].getInfo<CL_DEVICE_NAME>());
    }
    for (i = 0; i < gpuDevices.size(); i++)
    {   
        devices -> push_back(DeviceContainer(gpuDevices[i], *context));
        clExt = gpuDevices[i].getInfo<CL_DEVICE_EXTENSIONS>();
        doublePrecision -> push_back((clExt.find("cl_khr_fp64") != std::string::npos));
        deviceTypes -> push_back("GPU");
        deviceNames -> push_back(gpuDevices[i].getInfo<CL_DEVICE_NAME>());
    }
 
}

SpatialSEIR::PlatformContainer::~PlatformContainer()
{
    delete[] devices;
    delete[] deviceTypes;
    delete[] deviceNames;
    delete[] doublePrecision;
    delete context;
    delete platform;
}


SpatialSEIR::OCLProvider::OCLProvider()
{
    std::cout << "Setting up OpenCL Interface\n";
    try
    {
        // Allocate space for platforms and current config
        platforms = new std::vector<PlatformContainer>();
        currentPlatform = new cl::Platform*;
        currentContext = new cl::Context*;
        currentDevice = new DeviceContainer*;

        R_star_args = new FC_R_Star_KernelData(); 
        R_star_args -> totalWorkUnits = -1;

        // Allocate space for kernels
        test_kernel = new cl::Kernel();
        R_Star_kernel = new cl::Kernel();

        // Build platforms, devices, contexts
        cl_uint i;
        std::vector<cl::Platform> pformVec;
        cl::Platform::get(&pformVec);
        for (i = 0; i < pformVec.size(); i++) 
        {
            platforms -> push_back(PlatformContainer(pformVec[i]));
        }
        // Flag for existence of current<item>s
        isSetup = new int; *isSetup = 0;

    }
    catch(cl::Error e)
    {
        cout << "Problem getting platforms:" << endl;
        cout << e.what() << ": Error Code " << e.err() << endl;
        return;
    }

    // Create Kernels
    setDevice(0,0);
    *test_kernel = buildProgramForKernel("test_kernel.cl", (**currentDevice));
    *R_Star_kernel = buildProgramForKernel("R_Star_FC.cl", (**currentDevice));


    test();
}

SpatialSEIR::OCLProvider::~OCLProvider()
{
    delete currentPlatform;
    delete currentContext;
    delete currentDevice;
    delete platforms;
    delete R_star_args;
    delete test_kernel;
    delete R_Star_kernel;
    delete isSetup;
}

void SpatialSEIR::OCLProvider::setDevice(int platformId, int deviceId)
{
    *currentPlatform = (*platforms)[platformId].platform;
    *currentContext = (*platforms)[platformId].context;
    *currentDevice = &((*((*platforms)[platformId].devices))[deviceId]);
    *isSetup = 1;
}

cl::Kernel SpatialSEIR::OCLProvider::buildProgramForKernel(std::string kernelFile, DeviceContainer device)
{
    int err = 1;
    std::vector<cl::Device> devices; devices.push_back(*device.device);
    std::string log;
    // LKD is set at compile time, intall directory of OpenCL kernels. 
    std::string LKD(LSS_KERNEL_DIRECTORY);
    LKD = LKD.append(kernelFile);
    const char* progName = LKD.c_str();

    std::ifstream programFile(progName);
    std::string programString(std::istreambuf_iterator<char>(programFile), 
                             (std::istreambuf_iterator<char>()));
    cl::Program::Sources source(1, std::make_pair(programString.c_str(), 
                                programString.length() + 1));

    cl::Program program(**currentContext, source);
    std::vector<cl::Kernel> kernels;
    try
    {
        err = program.build(devices);
        log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
        program.createKernels(&kernels);
    }
    catch(cl::Error e)
    {
        std::cout << "CL Error in: " << e.what()<< "\n";
        std::cout << "CL Error: " << e.err()<< "\n";
        log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
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

    return(kernels[0]);
}



