#ifndef __CL_ENABLE_EXCEPTIONS
#define __CL_ENABLE_EXCEPTIONS
#endif

#include <CL/cl.hpp>
#include <OCLProvider.hpp>

namespace SpatialSEIR
{
    double OCLProvider::FC_R_Star_Part1(int startLoc,
                                 int startTime,
                                 int nTpts,
                                 int* R_star,
                                 int* S_star,
                                 int* R,
                                 int* I,
                                 double* p_rs,
                                 double p_ir)
    {
        cl::Device device = (cpuQueue -> getInfo<CL_QUEUE_DEVICE>());
        int remainingTpts = (nTpts - startTime);
        int i;
        void* mem_map;


        // Figure out a good way to partition the data 
        // and set up the index space. 
        //
        // Input:
        // 1. Integers (4 bytes)
        //    S_star (T)
        //    R_star (T)
        //    R      (T)
        //    I      (T)
        // 2. Doubles (8 bytes)
        //    p_rs   (T)
        size_t localMemPerCore = device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();
        int localSizeMultiple = (R_Star_p1_kernel -> getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device));
        int maxWorkUnits = (R_Star_p1_kernel -> getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device));
        int maxWorkUnitsCompile = (R_Star_p1_kernel -> getWorkGroupInfo<CL_KERNEL_COMPILE_WORK_GROUP_SIZE>(device))[0];
        int maxLocalSize = localMemPerCore/(4*4 + 8);
        int numWorkUnits = (maxLocalSize/localSizeMultiple)*localSizeMultiple; 
        int numWorkGroups = (remainingTpts/numWorkUnits); 
            numWorkGroups += (numWorkGroups*numWorkUnits < remainingTpts);
        int globalSize = numWorkGroups*numWorkUnits; 

        double* output = new double[numWorkGroups]();
        size_t buffSize = remainingTpts*sizeof(int);
        size_t localBuffSize = numWorkUnits*sizeof(int);

        std::cout << "localMemPerCore: " << localMemPerCore << "\n";
        std::cout << "localSizeMultiple: " << localSizeMultiple << "\n";
        std::cout << "maxWorkUnits: " << maxWorkUnits << "\n";
        std::cout << "maxWorkUnitsCompile: " << maxWorkUnitsCompile << "\n";
        std::cout << "maxLocalSize: " << maxLocalSize << "\n";
        std::cout << "numWorkUnits: " << numWorkUnits << "\n";
        std::cout << "numWorkGroups: " << numWorkGroups << "\n";
        std::cout << "globalSize: " << globalSize << "\n";


        cl::Buffer RstarBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_COPY_HOST_PTR, buffSize, R_star);
        cl::Buffer SstarBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_COPY_HOST_PTR, buffSize, S_star);
        cl::Buffer RBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_COPY_HOST_PTR, buffSize, R);
        cl::Buffer IBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_COPY_HOST_PTR, buffSize, I);
        cl::Buffer p_rsBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_COPY_HOST_PTR, remainingTpts*sizeof(double), p_rs);
        cl::Buffer outBuffer(*context, CL_MEM_READ_WRITE | 
            CL_MEM_COPY_HOST_PTR, numWorkGroups*sizeof(double), output);

        int err;  

        err = R_Star_p1_kernel -> setArg(0, RstarBuffer);
        err |= R_Star_p1_kernel -> setArg(1, SstarBuffer);
        err |= R_Star_p1_kernel -> setArg(2, RBuffer);
        err |= R_Star_p1_kernel -> setArg(3, IBuffer);
        err |= R_Star_p1_kernel -> setArg(4, p_rsBuffer);
        err |= R_Star_p1_kernel -> setArg(5, p_ir);
        err |= R_Star_p1_kernel -> setArg(6, remainingTpts);
        err |= R_Star_p1_kernel -> setArg(7, outBuffer);
        // Local Declarations
        err |= R_Star_p1_kernel -> setArg(8, localBuffSize, NULL); //R_star
        err |= R_Star_p1_kernel -> setArg(9, localBuffSize, NULL); //S_star
        err |= R_Star_p1_kernel -> setArg(10, localBuffSize, NULL); //R
        err |= R_Star_p1_kernel -> setArg(11, localBuffSize, NULL); //I
        err |= R_Star_p1_kernel -> setArg(12, numWorkUnits*sizeof(double), NULL); //p_rs

        if (err < 0)
        {
            std::cerr << "Couldn't set kernel args.\n";
            throw(err);
        }
        // Optimize this to use subbuffers so that we only write the data 
        // once per sampleR_star event
        try
        {
            cpuQueue -> enqueueNDRangeKernel(*R_Star_p1_kernel,
                                             0,
                                             globalSize,
                                             numWorkUnits,
                                             NULL,
                                             NULL
                                             );
        }
        catch(cl::Error e)
        {
            std::cerr << "Error enqueueing kernel: " << e.what() << "\n";
            std::cerr << "Error: " << e.err() << "\n";
            throw(-1);
        }

        mem_map = cpuQueue -> enqueueMapBuffer(outBuffer, CL_TRUE, CL_MAP_READ, 0, numWorkGroups*sizeof(double));
        memcpy(output, mem_map, numWorkGroups*sizeof(double));

        double outSum = 0.0;
        for (i = 0; i < numWorkGroups; i++)
        {
            outSum += output[i];
        }
        delete[] output;
        return(outSum);

    }

}
