#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void p_se_kernel1(int nLoc,
                           int nTpt,
                          __global int* I,
                          __global int* N,
                          __global double* eta,
                          __local int* I_loc,
                          __local int* N_loc,
                          __local double* eta_loc)
{
    size_t globalId = get_global_id(0);
    size_t localId = get_local_id(0);
    size_t localSize = get_local_size(0);
    size_t totalSize = nLoc*nTpt;

    if (globalId < totalSize)
    {
        I_loc[localId] = I[globalId];
        N_loc[localId] = N[globalId];
        eta_loc[localId] = eta[globalId];

        eta_loc[localId] = (exp(eta_loc[localId])*I_loc[localId])/(N_loc[localId]);
        
        eta[globalId] = eta_loc[localId];
    }
    else
    {
        eta_loc[localId] = 0.0;
    }

}


