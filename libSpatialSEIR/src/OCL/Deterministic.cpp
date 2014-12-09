#ifndef __CL_ENABLE_EXCEPTIONS
#define __CL_ENABLE_EXCEPTIONS
#endif

#ifdef ENABLE_OPEN_CL

#include <CL/cl.hpp>
#include <cmath>
#ifdef LSS_USE_BLAS
	#include <cblas.h>
#endif
#include <OCLProvider.hpp>
#include <ModelContext.hpp>
#include <CovariateMatrix.hpp>
#include <CompartmentalModelMatrix.hpp>
#include <DistanceMatrix.hpp>
#include <IOProvider.hpp>

#include <unistd.h>
namespace SpatialSEIR
{
    void OCLProvider::calculateP_SE(ModelContext* ctx)
    {
        cl::Context* context = *currentContext;
        cl::Device device = **((*currentDevice) -> device);

        cl::Event part1Finished, part2Finished, part3Finished;
        std::vector<cl::Event> waitList;
        waitList.push_back(part2Finished);

        ctx -> X -> calculate_eta_CPU(ctx -> eta, ctx -> beta);
        int nLoc = *(ctx -> S -> ncol);
        int nTpt = *(ctx -> S -> nrow); 

        size_t intBuffSize = nLoc*nTpt*sizeof(int);
        size_t doubleBuffSize = nLoc*nTpt*sizeof(double);
        cl::Buffer IBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_COPY_HOST_PTR, intBuffSize, (ctx -> I -> data));
        cl::Buffer NBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_COPY_HOST_PTR, intBuffSize, (ctx -> N));

        cl::Buffer offsetBuffer(*context, CL_MEM_READ_ONLY |
                CL_MEM_COPY_HOST_PTR, nTpt*sizeof(double), ctx -> offset);
        cl::Buffer etaBuffer(*context, CL_MEM_READ_WRITE | 
            CL_MEM_COPY_HOST_PTR, doubleBuffSize, (ctx -> eta));
        cl::Buffer p_seBuffer(*context, CL_MEM_READ_WRITE | 
            CL_MEM_COPY_HOST_PTR, doubleBuffSize, (ctx -> p_se));
        std::vector<cl::Buffer> distMatBuffers;
        int i;
        for (i = 0; i < (int) (ctx -> scaledDistMatrices -> size()); i++) 
        {
            distMatBuffers.push_back(
                cl::Buffer(*context, CL_MEM_WRITE_ONLY | 
                        CL_MEM_COPY_HOST_PTR, nLoc*nLoc*sizeof(double) , 
                        (*(ctx -> scaledDistMatrices))[i] -> data)
                );
        }


        // Kernel 1
        // Input:
        // 1. Integers (4 bytes)
        //    I (TxP)
        //    N (TxP)
        // 2. Doubles (8 bytes)
        //    eta (not exponentiated) (TxP)
        int workGroupSize;
        int numWorkGroups;
        int totalWorkUnits;
        int globalSize;

        if (R_star_args -> totalWorkUnits == -1)
        {
            // Not populated, need to calculate kernel parameters. 
            totalWorkUnits = nLoc*nTpt;
            size_t localMemPerCore = device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();
            int deviceMaxSize = (device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());
            int localSizeMultiple = (p_se_kernel1 -> getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device));
            int reportedMaxSize = (p_se_kernel1 -> getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device));
            int maxLocalSize = localMemPerCore/(2*4 + 1*8);
            maxLocalSize = std::min(std::min(maxLocalSize, deviceMaxSize), reportedMaxSize);
            int i=-1;
            workGroupSize = 0;
            while(workGroupSize < maxLocalSize && workGroupSize < totalWorkUnits)
            {
                i++;
                workGroupSize = pow(2,i)*localSizeMultiple;
            }
            if (workGroupSize >= maxLocalSize)
            {
                workGroupSize = pow(2,i-1)*localSizeMultiple;
            }

            numWorkGroups = (totalWorkUnits/workGroupSize); 
                numWorkGroups += (numWorkGroups*workGroupSize < totalWorkUnits);
            globalSize = numWorkGroups*workGroupSize;

            (p_se_args -> workGroupSize) = workGroupSize; 
            (p_se_args -> numWorkGroups) = numWorkGroups; 
            (p_se_args -> totalWorkUnits) = totalWorkUnits; 
            (p_se_args -> globalSize) = globalSize; 

            /*
            lssCout << "Total Work Units: " << totalWorkUnits << "\n";
            lssCout << "Local Mem Per Core: " << localMemPerCore << "\n";
            lssCout << "Local Size Multiple: " << localSizeMultiple << "\n";
            lssCout << "Reported Maximum Size: " << deviceMaxSize << "\n";
            lssCout << "Max Local Size: " << maxLocalSize << "\n";
            lssCout << "Work Group Size: " << workGroupSize << "\n";
            lssCout << "Global Size: " << globalSize << "\n";
            */
        }
        else
        {
            workGroupSize = (p_se_args -> workGroupSize); 
            numWorkGroups = (p_se_args -> numWorkGroups); 
            totalWorkUnits = (p_se_args -> totalWorkUnits); 
            globalSize = (p_se_args -> globalSize); 
        }

        int err;

        try
        {
            err = p_se_kernel1 -> setArg(0, nLoc);
            err |= p_se_kernel1 -> setArg(1, nTpt);
            err |= p_se_kernel1 -> setArg(2, IBuffer);
            err |= p_se_kernel1 -> setArg(3, NBuffer);
            err |= p_se_kernel1 -> setArg(4, etaBuffer);
            err |= p_se_kernel1 -> setArg(5, workGroupSize*sizeof(int), NULL); //I
            err |= p_se_kernel1 -> setArg(6, workGroupSize*sizeof(int), NULL); //N
            err |= p_se_kernel1 -> setArg(7, workGroupSize*sizeof(double), NULL); //eta
            if (err < 0)
            {
                std::cerr << "Couldn't set kernel args.\n";
                throw(err);
            }
        }
        catch(cl::Error e)
        {
            std::cerr << "Error setting kernel arguments: " << e.what() << "\n";
            std::cerr << "Error: " << e.err() << "\n";
            throw(-1);
        }
        try
        {
            ((*currentDevice) -> commandQueue) -> enqueueNDRangeKernel(*p_se_kernel1,
                                                                     cl::NDRange(0),
                                                                     cl::NDRange(globalSize),
                                                                     cl::NDRange(workGroupSize),
                                                                     NULL,
                                                                     &part1Finished
                                                                     );
        }
        catch(cl::Error e)
        {
            std::cerr << "Error enqueueing kernel: " << e.what() << "\n";
            std::cerr << "Error: " << e.err() << "\n";
            throw(-1);
        }
    
        // Kernel 2
        // Input:
        // 1. Doubles (8 bytes)
        //    p_se_components (TxP) (already on device)
        //    scaled distance matrix (PxP)
        //    p_se for output matrix (TxP)



        // Translation from pure c++ process:
        // 1. Exponentiate eta
        // 2. eta -> p_seComponents
        // 3. Matrix mult
        //    (A) p_se as output. 
        //    (B) p_se_components as first matrix
        //    (C) scaled distance matrix as second matrix
        // 4. Final assembly


        /*
         * clblasOrder order
         * clblasTranspose transA
         * clblasTranspose transB
         * size_t M,
         * size_t N,
         * size_t K,
         * cl_double alpha,
         * const cl_mem A,
         * size_t offA,
         * size_t lda,
         * const cl_mem B,
         * size_t offB,
         * size_t ldb,
         * cl_double beta,
         * cl_mem C,
         * size_t offC,
         * size_t ldc,
         * cl_uint numCommandQueues,
         * cl_command_queue* commandQueues,
         * cl_uint numEventsInWaitList,
         * const cl_event* eventWaitList,
         * cl_event* events
         */
        cl_uint numCommandQueues = 1;
        
        try
        {
            int numMatrices = (ctx -> scaledDistMatrices -> size());
            for (i = 0; i < numMatrices; i++)
            {
               clblasStatus multErr = clblasDgemm(clblasColumnMajor,   // Order
                                               clblasNoTrans,        // TransB
                                               clblasNoTrans,        // TransA
                                               nTpt,                 // M
                                               nLoc,                 // N
                                               nLoc,                 // K
                                               (ctx -> rho)[i],      // alpha
                                               etaBuffer(),          // A
                                               0,                    // offA
                                               nTpt,                 // ldA
                                               (distMatBuffers[i])(),// B
                                               0,                    // offB
                                               nLoc,                 // ldB
                                               (i!=0),               // Beta
                                               p_seBuffer(),         // C
                                               0,                    // offC
                                               nTpt,                 // ldC
                                               numCommandQueues,     // numCommandQueues
                                               &((*(**currentDevice).commandQueue)()), // commandQueues
                                               1,                    // numEventsInWaitList
                                               &(part1Finished()),   // eventWaitList
                                               &((waitList[0])()));  // events 
               
                if (multErr != CL_SUCCESS)
                { 
                    lssCout << "clBLAS Error Encountered: " << multErr << "\n";
                    throw(-1);
                }
            }
            
        }
        catch(cl::Error e)
        {
            throw(-1);
        }

        // Kernel 3
        // Input:
        // 1. Doubles (8 bytes)
        //    p_se_components (TxP) (already on device)
        //    p_se for output matrix (TxP) (already on device)
        //
        //

        // Part 3 kernel
        try
        {
            err = p_se_kernel2 -> setArg(0, nLoc);
            err |= p_se_kernel2 -> setArg(1, nTpt);
            err |= p_se_kernel2 -> setArg(2, etaBuffer);
            err |= p_se_kernel2 -> setArg(3, offsetBuffer);
            err |= p_se_kernel2 -> setArg(4, p_seBuffer);
            err |= p_se_kernel2 -> setArg(5, workGroupSize*sizeof(double), NULL); // p_seComponents loc
            err |= p_se_kernel2 -> setArg(6, workGroupSize*sizeof(double), NULL); // p_se loc
            if (err < 0)
            {
                lssCout << "Error setting kernel args \n";
            }
        }
        catch(cl::Error e)
        {
            std::cerr << "Error setting kernel arguments: " << e.what() << "\n";
            std::cerr << "Error: " << e.err() << "\n";
            throw(-1);
        }

        try
        {
            ((*currentDevice) -> commandQueue) -> enqueueNDRangeKernel(*p_se_kernel2,
                                                                     cl::NDRange(0),
                                                                     cl::NDRange(globalSize),
                                                                     cl::NDRange(workGroupSize),
                                                                     &waitList,
                                                                     &part3Finished
                                                                     );
            clWaitForEvents(1, &(part3Finished()));

        }
        catch(cl::Error e)
        {
            std::cerr << "Error enqueueing kernel: " << e.what() << "\n";
            std::cerr << "Error: " << e.err() << "\n";
            throw(-1);
        }
        void* p_seMap = ((*currentDevice) -> commandQueue) -> enqueueMapBuffer(
                p_seBuffer, CL_TRUE, CL_MAP_READ, 0, totalWorkUnits*sizeof(double));
        memcpy(ctx -> p_se, p_seMap, totalWorkUnits*sizeof(double));
        ((*currentDevice) -> commandQueue) -> enqueueUnmapMemObject(p_seBuffer, p_seMap);
        return;
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

        cl::Buffer bufferA(**currentContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 100*sizeof(float), A);
        cl::Buffer bufferB(**currentContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 100*sizeof(float), B);
        cl::Buffer bufferOut(**currentContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 100*sizeof(float), out);

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
            err = (*currentDevice) -> commandQueue -> enqueueNDRangeKernel(*test_kernel, 
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

        mappedMemory = ((*currentDevice) -> commandQueue -> enqueueMapBuffer(bufferOut, CL_TRUE, CL_MAP_READ, 0, 100*sizeof(float)));


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
        ((*currentDevice) -> commandQueue) -> enqueueUnmapMemObject(bufferOut, mappedMemory);

        delete[] A;
        delete[] B;
        delete[] out;
        delete[] cpuOut;
        return(0); 
    }


}

#endif

