#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void p_se_kernel2(int nLoc,
                           int nTpt,
                           double rho,
                          __global double* p_se_components,
                          __global double* p_se,
                          __local double* p_se_components_loc,
                          __local double* p_se_loc)

{
    size_t globalId = get_global_id(0);
    size_t localId = get_local_id(0);
    size_t localSize = get_local_size(0);
    size_t totalSize = nLoc*nTpt;

    if (globalId < totalSize)
    {
        p_se_components_loc[localId] = p_se_components[globalId];
        p_se_loc[localId] = p_se[globalId]; 
        p_se[globalId] = 1 - exp(-p_se_components_loc[localId] - rho*p_se_loc[localId]);  
    }
}


