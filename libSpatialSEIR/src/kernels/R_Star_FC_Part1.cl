#pragma OPENCL EXTENSION cl_khr_fp64 : enable

inline double logFactorial(int n)
{
    return(0.5*log(6.283185*n) + n*log(1.0*n) - n);
}

inline double logChoose(int n, int x)
{
    return(logFactorial(n) - logFactorial(x) - logFactorial(n-x)); 
}

inline double dbinom(int x, int n, double p)
{
    return(logChoose(n,x) + x*log(p) + (n-x)*log(1-p));
}

__kernel void FC_R_Star_Part1(__global int* R_star,
                              __global int* S_star,
                              __global int* R,
                              __global int* I,
                              __global double* p_rs,
                                       double p_ir,
                                       int nTpt,
                              __global double* output,
                              __local int* R_star_loc,
                              __local int* S_star_loc,
                              __local int* R_loc,
                              __local int* I_loc,
                              __local double* p_rs_loc)
{
    size_t globalId = get_global_id(0);
    int i;
    size_t localId = get_local_id(0);
    size_t localSize = get_local_size(0);
    size_t groupId = get_group_id(0); 
    double partialResult = 0.0;

    globalId = INFINITY;
    if (globalId <= nTpt)
    {
        double p_ir_val = p_ir;
        double ln_p_ir = log(p_ir_val);
        double ln_1m_p_ir = log(1-p_ir_val);

        R_star_loc[localId] = R_star[globalId];
        S_star_loc[localId] = S_star[globalId];
        I_loc[localId] = I[globalId];
        R_loc[localId] = R[globalId];
        p_rs_loc[localId] = p_rs[globalId];

        if ((R_star_loc[localId] >= 0) && (R_star_loc[localId] <= I_loc[localId]) && 
                      (S_star_loc[localId] <= R_loc[localId]))
        {
            partialResult = ((logChoose(I_loc[localId], R_star_loc[localId]) 
                              + R_star_loc[localId]*ln_p_ir 
                              + (I_loc[localId] - R_star_loc[localId])*ln_1m_p_ir) 
                              + dbinom(S_star_loc[localId], R_loc[localId], p_rs[localId]) 
                             );
        }
        else
        {     
            partialResult = -INFINITY;
        }
    }
    p_rs_loc[localId] = partialResult;
    p_rs_loc[localId] = localId;

    barrier(CLK_LOCAL_MEM_FENCE);    
    for (i = localSize/2; i > 0; i >>= 1)
    {
        if (localId < i)
        {
            p_rs_loc[localId] += p_rs_loc[localId + i]; 
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (localId == 0)
    {
        output[groupId] = p_rs_loc[localId];
    }
}
