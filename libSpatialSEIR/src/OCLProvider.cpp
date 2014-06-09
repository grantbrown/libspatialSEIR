#define __CL_ENABLE_EXCEPTIONS

#include <CL/cl.hpp>
#include <OCLProvider.hpp>
#include<iostream>
#include<fstream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>


SpatialSEIR::DeviceContainer::DeviceContainer(cl::Device *inDevice, cl::Context *inContext)
{
    device = new cl::Device*;
    *device = inDevice; 
    commandQueue = new cl::CommandQueue(*inContext, *inDevice); 
}
SpatialSEIR::DeviceContainer::~DeviceContainer()
{
    delete commandQueue;
    delete device;
}

SpatialSEIR::PlatformContainer::PlatformContainer(cl::Platform *inPlatform)
{
    //Store platform
    platform = new cl::Platform*; 
    *platform = inPlatform;

    //Store all platform devices to create context
    std::vector<cl::Device>* allDevices = new std::vector<cl::Device>();
    (*platform) -> getDevices(CL_DEVICE_TYPE_ALL, allDevices);
    context = new cl::Context(*allDevices);

    // Process and store devices by type
    devices = new std::vector<DeviceContainer*>();
    deviceTypes = new std::vector<std::string>();
    deviceNames = new std::vector<std::string>();
    doublePrecision = new std::vector<cl_uint>();

    std::vector<cl::Device> *cpuDevices = new std::vector<cl::Device>; 
    std::vector<cl::Device> *gpuDevices = new std::vector<cl::Device>; 

    try
    {
        (*platform) -> getDevices(CL_DEVICE_TYPE_CPU, cpuDevices);
    }
    catch(cl::Error e)
    {
       if (e.err() == -1) 
       {
           std::cout << "...no OpenCL CPU devices found\n";
       }
       else
       {
           std::cerr << "Error querying CPU devices.\n";
           std::cerr << e.what() << "\n";
           throw(e.err());
       }

    }

    try
    {
        (*platform) -> getDevices(CL_DEVICE_TYPE_GPU, gpuDevices);
    }
    catch(cl::Error e)
    {
       if (e.err() == -1) 
       {
           std::cout << "...no OpenCL GPU devices found\n";
       }
       else
       {
           std::cerr << "Error querying GPU devices.\n";
           std::cerr << e.what() << "\n";
           throw(e.err());
       }
    }


    std::string clExt;
    cl_uint i;
    DeviceContainer* newDevice;
    for (i = 0; i < cpuDevices -> size(); i++)
    {  
        std::cout << "Adding CPU Device: ";
        newDevice = new DeviceContainer(&(*cpuDevices)[i], context);
        devices -> push_back(newDevice);
        clExt = (*(newDevice -> device)) -> getInfo<CL_DEVICE_EXTENSIONS>();
        doublePrecision -> push_back((clExt.find("cl_khr_fp64") != std::string::npos));
        deviceTypes -> push_back("CPU");
        deviceNames -> push_back((*(newDevice -> device)) -> getInfo<CL_DEVICE_NAME>());
        std::cout << (*deviceNames)[deviceNames -> size() - 1] << "\n";

    }
    for (i = 0; i < gpuDevices -> size(); i++)
    {   
        std::cout << "Adding CPU Device: ";
        newDevice = new DeviceContainer(&((*gpuDevices)[i]), context);
        devices -> push_back(newDevice);
        clExt = (*(newDevice -> device)) -> getInfo<CL_DEVICE_EXTENSIONS>();
        doublePrecision -> push_back((clExt.find("cl_khr_fp64") != std::string::npos));
        deviceTypes -> push_back("GPU");
        deviceNames -> push_back((*(newDevice -> device)) -> getInfo<CL_DEVICE_NAME>());
        std::cout << (*deviceNames)[deviceNames -> size() - 1] << "\n";
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
        platforms = new std::vector<PlatformContainer*>();
        currentPlatform = new cl::Platform*;
        currentContext = new cl::Context*;
        currentDevice = new DeviceContainer*;

        R_star_args = new FC_R_Star_KernelData(); 
        R_star_args -> totalWorkUnits = -1;

        // Allocate space for kernels
        test_kernel = new cl::Kernel();
        R_Star_kernel = new cl::Kernel();
        p_se_kernel1 = new cl::Kernel();
        p_se_kernel2 = new cl::Kernel();


        // Build platforms, devices, contexts
        cl_uint i;
        std::vector<cl::Platform> *pformVec = new std::vector<cl::Platform>;
        cl::Platform::get(pformVec);

        PlatformContainer* newPlatform;
        for (i = 0; i < pformVec -> size(); i++) 
        {
            newPlatform = new PlatformContainer((&(*pformVec)[i]));
            platforms -> push_back(newPlatform);
        }


        // Initialize clBLAS library
        
        clblasStatus err = clblasSetup();
        if (err != CL_SUCCESS)
        {
            std::cout << "Error setting up clBLAS library: " << err << "\n";
            throw(-1);
        }

        // Flag for existence of current<item>s
        isSetup = new int; *isSetup = 0;


    }
    catch(cl::Error e)
    {
        cout << "Problem getting platforms:" << endl;
        cout << e.what() << ": Error Code " << e.err() << endl;
        throw(-1);
    }

    // Create Kernels
    // Dummy code to pick device 0,0
    setDevice(0,0);

    *test_kernel = buildProgramForKernel("test_kernel.cl", (*currentDevice));

    *R_Star_kernel = buildProgramForKernel("R_Star_FC.cl", (*currentDevice));

    *p_se_kernel1 = buildProgramForKernel("p_se_kernel1.cl", (*currentDevice));

    *p_se_kernel2 = buildProgramForKernel("p_se_kernel2.cl", (*currentDevice));

    test();
}

SpatialSEIR::OCLProvider::~OCLProvider()
{
    clblasTeardown();
    delete currentPlatform;
    delete currentContext;
    delete currentDevice;
    delete platforms;
    delete R_star_args;
    delete test_kernel;
    delete R_Star_kernel;
    delete p_se_kernel1;
    delete p_se_kernel2;
    delete isSetup;
}

void SpatialSEIR::OCLProvider::setDevice(int platformId, int deviceId)
{
    *currentPlatform = (*((*platforms)[platformId] -> platform));
    *currentContext = (*platforms)[platformId] -> context;
    *currentDevice = ((*((*platforms)[platformId] -> devices))[deviceId]);
    *isSetup = 1;
}

cl::Kernel SpatialSEIR::OCLProvider::buildProgramForKernel(std::string kernelFile, DeviceContainer* device)
{
    int err = 1;
    std::vector<cl::Device> devices; devices.push_back(**(device -> device));
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



