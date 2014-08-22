#define __CL_ENABLE_EXCEPTIONS

#include<math.h>
#include<cstring>
#include<vector>
#include <ModelContext.hpp>
#include <OCLProvider.hpp>
#include <IOProvider.hpp>

#ifdef ENABLE_OPEN_CL
#include <CL/cl.hpp>

SpatialSEIR::DeviceContainer::DeviceContainer(cl::Device *inDevice, cl::Context *inContext)
{
    device = new cl::Device*;
    *device = inDevice; 
    commandQueue = new cl::CommandQueue(*inContext, *inDevice); 
}
SpatialSEIR::DeviceContainer::~DeviceContainer()
{
    lssCout << "Deleting DeviceContainer.\n";
    delete commandQueue;
    lssCout << "1\n";
    delete device;
    lssCout << "2\n";
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
           lssCout << "...no OpenCL CPU devices found\n";
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
           lssCout << "...no OpenCL GPU devices found\n";
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
        lssCout << "Adding CPU Device: ";
        newDevice = new DeviceContainer(&(*cpuDevices)[i], context);
        devices -> push_back(newDevice);
        clExt = (*(newDevice -> device)) -> getInfo<CL_DEVICE_EXTENSIONS>();
        doublePrecision -> push_back((clExt.find("cl_khr_fp64") != std::string::npos));
        deviceTypes -> push_back("CPU");
        deviceNames -> push_back((*(newDevice -> device)) -> getInfo<CL_DEVICE_NAME>());
        lssCout << (*deviceNames)[deviceNames -> size() - 1] << "\n";

    }
    for (i = 0; i < gpuDevices -> size(); i++)
    {   
        lssCout << "Adding GPU Device: ";
        newDevice = new DeviceContainer(&((*gpuDevices)[i]), context);
        devices -> push_back(newDevice);
        clExt = (*(newDevice -> device)) -> getInfo<CL_DEVICE_EXTENSIONS>();
        doublePrecision -> push_back((clExt.find("cl_khr_fp64") != std::string::npos));
        deviceTypes -> push_back("GPU");
        deviceNames -> push_back((*(newDevice -> device)) -> getInfo<CL_DEVICE_NAME>());
        lssCout << (*deviceNames)[deviceNames -> size() - 1] << "\n";
    }
}

SpatialSEIR::PlatformContainer::~PlatformContainer()
{
    while (devices -> size() != 0){delete (devices -> back()); devices -> pop_back();}
    delete devices;
    delete deviceTypes;
    while (deviceNames -> size() != 0){deviceNames -> pop_back();}
    delete deviceNames;
    while (doublePrecision -> size() != 0){doublePrecision -> pop_back();}
    delete doublePrecision;
    delete context;
    delete platform;
}


SpatialSEIR::OCLProvider::OCLProvider()
{
    lssCout << "Setting up OpenCL Interface\n";
    try
    {
        // Allocate space for platforms and current config
        platforms = new std::vector<PlatformContainer*>();
        currentPlatform = new cl::Platform*;
        currentContext = new cl::Context*;
        currentDevice = new DeviceContainer*;

        R_star_args = new FC_R_Star_KernelData(); 
        R_star_args -> totalWorkUnits = -1;

        p_se_args = new P_SE_Calculation_KernelData();
        p_se_args -> totalWorkUnits = -1;

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
            lssCout << "Error setting up clBLAS library: " << err << "\n";
            throw(-1);
        }

        // Flag for existence of current<item>s
        isSetup = new int; *isSetup = 0;


    }
    catch(cl::Error e)
    {
        lssCout << "Problem getting platforms:\n";
        lssCout << e.what() << ": Error Code " << e.err() << "\n";
        throw(-1);
    }

    // Create Kernels
    // Dummy code to pick device 0,0
    setDevice(0,0);
    buildKernels();
    test();
}

void SpatialSEIR::OCLProvider::buildKernels()
{
    *test_kernel = buildProgramForKernel("test_kernel.cl", (*currentDevice));
    *R_Star_kernel = buildProgramForKernel("R_Star_FC.cl", (*currentDevice));
    *p_se_kernel1 = buildProgramForKernel("p_se_kernel1.cl", (*currentDevice));
    *p_se_kernel2 = buildProgramForKernel("p_se_kernel2.cl", (*currentDevice));
}

void SpatialSEIR::OCLProvider::printSummary()
{
    unsigned int i,j;
    if (platforms -> size() == 0)
    {
        lssCout << "No OpenCL Platforms Detected\n";
        return;
    }
    for (i = 0; i < (platforms -> size()); i++)
    {
        lssCout << "Platform " << i << ": "  << (*(((*platforms)[i]) -> platform)) -> getInfo<CL_PLATFORM_NAME>() << "\n";
        lssCout << "  Devices: \n";
        for (j = 0; j < (((*(platforms))[i]) -> devices) -> size(); j++)
        {
            lssCout << "  " << j << ". " <<  (*(((*(((*(platforms))[i]) -> devices))[j]) -> device)) -> getInfo<CL_DEVICE_NAME>() << "\n";
            lssCout << (j < 10 ? " " : (j < 100 ? "  " : (j < 1000 ? "   " : "   "))) << "    " << "Type: " << 
                ((*(((*(platforms))[i]) -> deviceTypes))[j])  << "\n";
            lssCout << (j < 10 ? " " : (j < 100 ? "  " : (j < 1000 ? "   " : "   "))) << "    " << "Supports Double Precision: " << 
                ((*(((*(platforms))[i]) -> doublePrecision))[j])  << "\n";
            lssCout << (j < 10 ? " " : (j < 100 ? "  " : (j < 1000 ? "   " : "   "))) << "    " << 
                "Preferred double vector width: " << (*(((*(((*(platforms))[i]) -> devices))[j]) -> device)) -> 
                getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE>() << "\n";
            lssCout << (j < 10 ? " " : (j < 100 ? "  " : (j < 1000 ? "   " : "   "))) << "    " << 
                "Local Memory: " << (*(((*(((*(platforms))[i]) -> devices))[j]) -> device)) -> getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << "\n";
            lssCout << (j < 10 ? " " : (j < 100 ? "  " : (j < 1000 ? "   " : "   "))) << "    " << 
                "Max Compute Units: " << (*(((*(((*(platforms))[i]) -> devices))[j]) -> device)) -> getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << "\n";
        }

    }
}

SpatialSEIR::OCLProvider::~OCLProvider()
{
    clblasTeardown();
    delete currentPlatform;
    delete currentContext;
    delete currentDevice;
    delete platforms;
    delete R_star_args;
    delete p_se_args;
    delete test_kernel;
    delete R_Star_kernel;
    delete p_se_kernel1;
    delete p_se_kernel2;
    delete isSetup;
}

void SpatialSEIR::OCLProvider::setDevice(int platformId, int deviceId)
{

    unsigned int pID, dID;
    if (platformId < 0 || deviceId < 0)
    {
        std::cerr << "Invalid Arguments\n";
    }
    pID = platformId;
    dID = deviceId;
    if ((*platforms).size() < pID)
    {
         std::cerr << "Invalid Platform.\n";
         throw(-1);      
    }
    if ((*((*platforms)[pID] -> devices)).size() <= dID)
    {
        std::cerr << "Invalid Device.\n";
        throw(-1);
    }

    // Make sure R_star local size is re-calculated
    R_star_args -> totalWorkUnits = -1;
    p_se_args -> totalWorkUnits = -1;

    *currentPlatform = (*((*platforms)[pID] -> platform));
    *currentContext = (*platforms)[pID] -> context;
    *currentDevice = ((*((*platforms)[pID] -> devices))[dID]);
    buildKernels();
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
        if (log.find("warning") != std::string::npos)
        {
            lssCout << "Warnings generated while building kernel.\n";
            lssCout << "CL_PROGRAM_BUILD_LOG: \n" << log << "\n";
        }
        program.createKernels(&kernels);
    }
    catch(cl::Error e)
    {
        lssCout << "CL Error in: " << e.what()<< "\n";
        lssCout << "CL Error: " << e.err()<< "\n";
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
#endif
#ifndef ENABLE_OPEN_CL

SpatialSEIR::OCLProvider::OCLProvider()
{
    // Do Nothing
}

void SpatialSEIR::OCLProvider::setDevice(int platformId, int deviceId)
{
    lssCout << "OpenCL not supported.\n";
    throw(-1);
}

void SpatialSEIR::OCLProvider::printSummary()
{
    lssCout << "OpenCL not supported.\n";
    throw(-1);
}

int SpatialSEIR::OCLProvider::test()
{
    lssCout << "OpenCL not supported.\n";
    throw(-1);
}

double SpatialSEIR::OCLProvider::FC_R_Star(int nLoc, int nTpts, int* S_star, 
                                    int* E_star, int* R_star, int* S,
                                    int* I, int* R, double* p_se,
                                    double* p_rs, double* p_ir)
{
    lssCout << "OpenCL not supported.\n";
    throw(-1);
}

void SpatialSEIR::OCLProvider::calculateP_SE(ModelContext* ctx)
{
    lssCout << "OpenCL not supported.\n";
    throw(-1);
}

SpatialSEIR::OCLProvider::~OCLProvider()
{
    // Do nothing
}







#endif


