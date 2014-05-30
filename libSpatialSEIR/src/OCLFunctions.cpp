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
        int remainingTpts = (nTpts - startTime);
        double* output = new double[remainingTpts]();

        cl::Buffer RstarBuffer(*context, CL_MEM_READ_WRITE | 
            CL_MEM_COPY_HOST_PTR, remainingTpts*sizeof(int), R_star);
        cl::Buffer SstarBuffer(*context, CL_MEM_READ_WRITE | 
            CL_MEM_COPY_HOST_PTR, remainingTpts*sizeof(int), S_star);
        cl::Buffer RBuffer(*context, CL_MEM_READ_WRITE | 
            CL_MEM_COPY_HOST_PTR, remainingTpts*sizeof(int), R);
        cl::Buffer IBuffer(*context, CL_MEM_READ_WRITE | 
            CL_MEM_COPY_HOST_PTR, remainingTpts*sizeof(int), I);
        cl::Buffer p_rsBuffer(*context, CL_MEM_READ_WRITE | 
            CL_MEM_COPY_HOST_PTR, remainingTpts*sizeof(double), p_rs);
        cl::Buffer outBuffer(*context, CL_MEM_READ_ONLY | 
            CL_MEM_COPY_HOST_PTR, remainingTpts*sizeof(double), output);

        int err;  

        err = R_Star_p1_kernel -> setArg(0, RstarBuffer);
        err |= R_Star_p1_kernel -> setArg(1, SstarBuffer);
        err |= R_Star_p1_kernel -> setArg(2, RBuffer);
        err |= R_Star_p1_kernel -> setArg(3, IBuffer);
        err |= R_Star_p1_kernel -> setArg(4, p_rsBuffer);
        err |= R_Star_p1_kernel -> setArg(5, p_ir);
        err |= R_Star_p1_kernel -> setArg(6, output);

        if (err < 0)
        {
            std::cerr << "Couldn't set kernel args.\n";
            throw(err);
        }

        // Optimize this to use subbuffers so that we only write the data 
        // once per sampleR_star event
        cpuQueue -> enqueueNDRangeKernel(*R_Star_p1_kernel,
                                         0,
                                         remainingTpts,
                                         1,
                                         NULL,
                                         NULL
                                         );

        // Do stuff with output

        delete[] output;

        return(0.0);

    }

}
