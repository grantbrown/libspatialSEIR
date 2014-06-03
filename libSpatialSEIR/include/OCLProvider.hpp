#ifndef OCL_PROVIDER_INC
#define OCL_PROVIDER_INC

#define __CL_ENABLE_EXCEPTIONS

#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#endif

#include <CL/cl.hpp>


namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class OCLProvider
    {
        public:  
            //Methods 
            OCLProvider();

            std::vector<cl::Kernel>* buildProgramForKernel(std::string kernelFile, 
                        std::vector<cl::Device> devices);
            int test();

            ~OCLProvider();

            //Attributes
            cl::Context *context;
            std::vector<cl::Platform> *platforms;
            std::vector<cl::Device> *platformDevices, 
                                    *allDevices, 
                                    *ctxDevices; 
            std::vector<std::string> *deviceNames;
            std::vector<std::vector<size_t> > *workSizes;
            std::vector<size_t> *globalMemSizes;
            std::vector<cl_ulong> *localMemSizes;
            std::vector<size_t> *numComputeUnits;
            std::vector<cl_uint> *doublePrecision;
            std::vector<cl::Program> *programs;

            // FC Methods
            double FC_R_Star_Part1(int nLoc, 
                                   int nTpts,
                                   int* R_star,
                                   int* S_star,
                                   int* R,
                                   int* I,
                                   double* p_rs,
                                   double p_ir);

            // Kernels
            cl::Kernel* test_kernel;
            cl::Kernel* R_Star_p1_kernel;

            // Queues
            cl::CommandQueue* cpuQueue;
            cl::CommandQueue* gpuQueue;
    };
}

#endif
