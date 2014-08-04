#ifndef OCL_PROVIDER_INC
#define OCL_PROVIDER_INC

#define __CL_ENABLE_EXCEPTIONS

#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<math.h>
#include<cstring>
#include<vector>
#endif

#ifdef ENABLE_OPEN_CL

#include <CL/cl.hpp>
#include <ModelContext.hpp>
#include <clBLAS.h>

namespace SpatialSEIR
{

    struct P_SE_Calculation_KernelData
    {
        int workGroupSize;
        int numWorkGroups;
        int totalWorkUnits;
        int globalSize;
    };

    struct FC_R_Star_KernelData
    {
        int workGroupSize;
        int numWorkGroups;
        int totalWorkUnits;
        int globalSize;
        double *outCache;
    };

    class DeviceContainer
    {
        public:
            DeviceContainer(cl::Device *inDevice, cl::Context *inContext);
            ~DeviceContainer();
            cl::Device** device;
            cl::CommandQueue * commandQueue;
    };

    class PlatformContainer
    {
        public:
            PlatformContainer(cl::Platform *inPlatform);
            ~PlatformContainer();

            cl::Platform **platform;
            cl::Context *context;
            std::vector<DeviceContainer*> *devices;
            std::vector<std::string> *deviceTypes; // CPU or GPU
            std::vector<std::string> *deviceNames;
            std::vector<cl_uint> *doublePrecision;
    };

    class OCLProvider
    {
        public:  
            //Methods 
            OCLProvider();
            void setDevice(int platformId, int deviceId);
            void buildKernels();
            void printSummary();

            ~OCLProvider();

            //Attributes
            cl::Platform **currentPlatform;
            cl::Context **currentContext;
            DeviceContainer **currentDevice;
            int *isSetup;

            std::vector<PlatformContainer*> *platforms;
            cl::Kernel buildProgramForKernel(std::string kernelFile, 
                    DeviceContainer *device);
            std::vector<cl::Program> *programs;
            FC_R_Star_KernelData* R_star_args;
            P_SE_Calculation_KernelData* p_se_args;

            // FC Methods
            int test(); 
            double FC_R_Star(int nLoc, 
                             int nTpts,
                             int* S_star,
                             int* E_star,
                             int* R_star,
                             int* S,
                             int* I,
                             int* R,
                             double* p_se,
                             double* p_rs,
                             double* p_ir);

            void calculateP_SE(ModelContext* ctx);

            // Kernels
            cl::Kernel* test_kernel;
            cl::Kernel* R_Star_kernel;
            cl::Kernel* p_se_kernel1;
            cl::Kernel* p_se_kernel2;


            // Queues
            cl::CommandQueue* cpuQueue;
            cl::CommandQueue* gpuQueue;
    };
}

#endif


#ifndef ENABLE_OPEN_CL

namespace SpatialSEIR
{
    class OCLProvider
    {
        public:  
            //Methods 
            OCLProvider();
            void setDevice(int platformId, int deviceId);
            void printSummary();

            ~OCLProvider();

            // FC Methods
            int test(); 
            double FC_R_Star(int nLoc, 
                             int nTpts,
                             int* S_star,
                             int* E_star,
                             int* R_star,
                             int* S,
                             int* I,
                             int* R,
                             double* p_se,
                             double* p_rs,
                             double* p_ir);

            void calculateP_SE(ModelContext* ctx);
    };
}


#endif

#endif
