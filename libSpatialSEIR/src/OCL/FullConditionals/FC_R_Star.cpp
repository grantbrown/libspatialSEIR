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
                                  double* p_ir)
    {
        cl::Context* context = *currentContext;
        cl::Device device = **((*currentDevice) -> device);
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
            //    p_ir   (T)
            //    p_se   (TxP)
            R_star_args -> totalWorkUnits = nLoc*nTpts;
            size_t localMemPerCore = device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();
            int deviceMaxSize = (device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());
            int localSizeMultiple = (R_Star_kernel -> getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device));
            int reportedMaxSize = (R_Star_kernel -> getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device));


            int maxLocalSize = localMemPerCore/(4*6 + 3*8);
            maxLocalSize = std::min(std::min(maxLocalSize, deviceMaxSize), reportedMaxSize);
            i=-1;
            int workGroupSize = 0;
            while(workGroupSize <= maxLocalSize && workGroupSize < (R_star_args -> totalWorkUnits))
            {
                i++;
                workGroupSize = pow(2,i)*localSizeMultiple;
            }
            if (workGroupSize >= maxLocalSize)
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
        cl::Buffer p_irBuffer(*context, CL_MEM_WRITE_ONLY | 
            CL_MEM_USE_HOST_PTR, nTpts*sizeof(double), p_ir);
        cl::Buffer outBuffer(*context, CL_MEM_READ_WRITE | 
            CL_MEM_COPY_HOST_PTR, (R_star_args -> numWorkGroups)*sizeof(double), (R_star_args -> outCache));

        int err;  

        err = R_Star_kernel -> setArg(0, nTpts);
        err |= R_Star_kernel -> setArg(1, nLoc);
        err |= R_Star_kernel -> setArg(2, SstarBuffer);
        err |= R_Star_kernel -> setArg(3, EstarBuffer);
        err |= R_Star_kernel -> setArg(4, RstarBuffer);

        err |= R_Star_kernel -> setArg(5, SBuffer);
        err |= R_Star_kernel -> setArg(6, IBuffer);
        err |= R_Star_kernel -> setArg(7, RBuffer);
        err |= R_Star_kernel -> setArg(8, p_seBuffer);
        err |= R_Star_kernel -> setArg(9, p_rsBuffer);
        err |= R_Star_kernel -> setArg(10, p_irBuffer);
        err |= R_Star_kernel -> setArg(11, outBuffer);
        // Local Declarations
        err |= R_Star_kernel -> setArg(12, localBuffSize, NULL); //S_star
        err |= R_Star_kernel -> setArg(13, localBuffSize, NULL); //E_star
        err |= R_Star_kernel -> setArg(14, localBuffSize, NULL); //R_star
        err |= R_Star_kernel -> setArg(15, localBuffSize, NULL); //S
        err |= R_Star_kernel -> setArg(16, localBuffSize, NULL); //I
        err |= R_Star_kernel -> setArg(17, localBuffSize, NULL); //R
        err |= R_Star_kernel -> setArg(18, (R_star_args -> workGroupSize)*sizeof(double), NULL); //p_se
        err |= R_Star_kernel -> setArg(19, (R_star_args -> workGroupSize)*sizeof(double), NULL); //p_rs
        err |= R_Star_kernel -> setArg(20, (R_star_args -> workGroupSize)*sizeof(double), NULL); //p_ir

        if (err < 0)
        {
            std::cerr << "Couldn't set kernel args.\n";
            throw(err);
        }
        // Optimize this to use subbuffers so that we only write the data 
        // once per sampleR_star event
        try
        {
            (*currentDevice) -> commandQueue -> enqueueNDRangeKernel(*R_Star_kernel,
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

        mem_map = ((*currentDevice) -> commandQueue) -> enqueueMapBuffer(outBuffer, CL_TRUE, CL_MAP_READ, 0, (R_star_args -> numWorkGroups)*sizeof(double));
        memcpy((R_star_args -> outCache), mem_map, (R_star_args -> numWorkGroups)*sizeof(double));
        ((*currentDevice) -> commandQueue) -> enqueueUnmapMemObject(outBuffer, mem_map);

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
}

#endif

