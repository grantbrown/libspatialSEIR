#ifndef __CL_ENABLE_EXCEPTIONS
#define __CL_ENABLE_EXCEPTIONS
#endif

#include <CL/cl.hpp>
#include <cmath>
#include <OCLProvider.hpp>

namespace SpatialSEIR
{
    double OCLProvider::FC_R_Star(int nLoc, 
                                  int nTpts,
                                  int* S_star,
                                  int* E_star,
                                  int* R_star,
                                  int* S,
                                  int* I,
                                  int* R,
                                  double* p_se,
                                  double* p_rs,
                                  double p_ir)
    {
        cl::Device device = (cpuQueue -> getInfo<CL_QUEUE_DEVICE>());
        int i;       
        if (R_star_args -> totalWorkUnits == -1)
        {
            // Not populated, need to calculate kernel parameters. 
            
            // Figure out a good way to partition the data 
            // and set up the index space. 
            //
            // Input:
            // 1. Integers (4 bytes)
            //    S_star (TxP)
            //    E_star (TxP)
            //    R_star (TxP)
            //    S      (TxP)
            //    R      (TxP)
            //    I      (TxP)
            // 2. Doubles (8 bytes)
            //    p_rs   (T), might use as (TxP) for simplicity. 
            //    p_se   (TxP)
            R_star_args -> totalWorkUnits = nLoc*nTpts;
            size_t localMemPerCore = device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();
            int localSizeMultiple = (R_Star_p1_kernel -> getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device));
            int maxLocalSize = localMemPerCore/(4*6 + 2*8);
            i=-1;
            int workGroupSize = 0;
            while(workGroupSize <= maxLocalSize && workGroupSize < (R_star_args -> totalWorkUnits))
            {
                i++;
                workGroupSize = pow(2,i)*localSizeMultiple;
            }
            if (workGroupSize < (R_star_args -> totalWorkUnits))
            {
                workGroupSize = pow(2,i-1)*localSizeMultiple;
            }
            R_star_args -> workGroupSize = workGroupSize;

            int numWorkGroups = ((R_star_args -> totalWorkUnits)/workGroupSize); 
                numWorkGroups += (numWorkGroups*workGroupSize < (R_star_args -> totalWorkUnits));
            int globalSize = numWorkGroups*workGroupSize;
            double* output = new double[numWorkGroups]();

            R_star_args -> globalSize = globalSize;
            R_star_args -> numWorkGroups = numWorkGroups;
            R_star_args -> outCache = output; 
        }


        void* mem_map;
        size_t buffSize = (R_star_args -> totalWorkUnits)*sizeof(int);
        size_t localBuffSize = (R_star_args -> workGroupSize*sizeof(int));

        cl::Buffer SstarBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_USE_HOST_PTR, buffSize, S_star);
        cl::Buffer EstarBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_USE_HOST_PTR, buffSize, E_star);
        cl::Buffer RstarBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_USE_HOST_PTR, buffSize, R_star);
        cl::Buffer SBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_USE_HOST_PTR, buffSize, S);
        cl::Buffer IBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_USE_HOST_PTR, buffSize, I);
        cl::Buffer RBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_USE_HOST_PTR, buffSize, R);
        cl::Buffer p_seBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_USE_HOST_PTR, (R_star_args -> totalWorkUnits)*sizeof(double), p_se);
        cl::Buffer p_rsBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_USE_HOST_PTR, nTpts*sizeof(double), p_rs);
        cl::Buffer outBuffer(*context, CL_MEM_READ_WRITE | 
            CL_MEM_COPY_HOST_PTR, (R_star_args -> numWorkGroups)*sizeof(double), (R_star_args -> outCache));

        int err;  

        err = R_Star_p1_kernel -> setArg(0, nTpts);
        err |= R_Star_p1_kernel -> setArg(1, nLoc);
        err |= R_Star_p1_kernel -> setArg(2, SstarBuffer);
        err |= R_Star_p1_kernel -> setArg(3, EstarBuffer);
        err |= R_Star_p1_kernel -> setArg(4, RstarBuffer);

        err |= R_Star_p1_kernel -> setArg(5, SBuffer);
        err |= R_Star_p1_kernel -> setArg(6, IBuffer);
        err |= R_Star_p1_kernel -> setArg(7, RBuffer);
        err |= R_Star_p1_kernel -> setArg(8, p_seBuffer);
        err |= R_Star_p1_kernel -> setArg(9, p_rsBuffer);
        err |= R_Star_p1_kernel -> setArg(10, p_ir);
        err |= R_Star_p1_kernel -> setArg(11, outBuffer);
        // Local Declarations
        err |= R_Star_p1_kernel -> setArg(12, localBuffSize, NULL); //S_star
        err |= R_Star_p1_kernel -> setArg(13, localBuffSize, NULL); //E_star
        err |= R_Star_p1_kernel -> setArg(14, localBuffSize, NULL); //R_star
        err |= R_Star_p1_kernel -> setArg(15, localBuffSize, NULL); //S
        err |= R_Star_p1_kernel -> setArg(16, localBuffSize, NULL); //I
        err |= R_Star_p1_kernel -> setArg(17, localBuffSize, NULL); //R
        err |= R_Star_p1_kernel -> setArg(18, (R_star_args -> workGroupSize)*sizeof(double), NULL); //p_se
        err |= R_Star_p1_kernel -> setArg(19, (R_star_args -> workGroupSize)*sizeof(double), NULL); //p_rs

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
                                             R_star_args -> globalSize,
                                             R_star_args -> workGroupSize,
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

        mem_map = cpuQueue -> enqueueMapBuffer(outBuffer, CL_TRUE, CL_MAP_READ, 0, (R_star_args -> numWorkGroups)*sizeof(double));
        memcpy((R_star_args -> outCache), mem_map, (R_star_args -> numWorkGroups)*sizeof(double));
        cpuQueue -> enqueueUnmapMemObject(outBuffer, mem_map);

        double outSum = 0.0;
        for (i = 0; i < (R_star_args -> numWorkGroups); i++)
        {
            outSum += (R_star_args -> outCache)[i];
        } 

        if (std::isnan(outSum))
        {
            outSum = -INFINITY;
        }

        return(outSum);
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

        try
        {
            err = cpuQueue -> enqueueNDRangeKernel(*test_kernel, 
                                           0,           // Global Offset
                                           100,         // Global Size
                                           1,        // Local Size
                                           NULL,        // Event Vector
                                           NULL         // Event Pointer
                                           );
        }
        catch(cl::Error e)
        {
            std::cerr << "Error enqueueing kernel: " << e.what() << "\n";
            std::cerr << "Error: " << e.err() << "\n";
            throw(-1);
        }

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

}
