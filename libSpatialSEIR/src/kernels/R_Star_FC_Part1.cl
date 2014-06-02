#pragma OPENCL EXTENSION cl_khr_fp64 : enable
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
    if (globalId <= nTpt)
    {


        size_t localId = get_local_id(0);
        size_t localSize = get_local_size(0);
        size_t groupId = get_group_id(0); 

        double p_ir_val = p_ir;
        double ln_p_ir = log(p_ir_val);
        double ln_1m_p_ir = log(1-p_ir_val);

        R_star_loc[localId] = R_star[globalId];
        S_star_loc[localId] = S_star[globalId];
        I_loc[localId] = I[globalId];
        R_loc[localId] = R[globalId];
        p_rs_loc[localId] = p_rs[globalId];


        if (R_star_loc[localId] < 0 || R_star_loc[localId] > I_loc[localId] || 
                      S_star_loc[localId] > R_loc[localId])
        {
            output[globalId] = -INFINITY;
            return;
        }
        
        output[globalId] = (ln_p_ir*R_star_loc[localId] + 
                     ln_1m_p_ir*(I_loc[localId] - R_star_loc[localId]) + 
                     log(p_rs_loc[localId])*S_star_loc[localId] + log(1-p_rs_val[localId])*(R_loc[localId]) + 
                     R_loc[localId]*log(1.0*R_loc[localId]) - 
                     (R_loc[localId]-S_star_loc[localId])*log(1.0*R_loc[localId]-S_star_loc[localId]) - S_star_loc[localId] +  
                     I_loc[localId]*log(1.0*I_loc[localId]) - R_star_loc[localId]*log(1.0*R_star_loc[localId]) - 
                     (I_loc[localId] - R_star_loc[localId])*log(1.0*I_loc[localId] - R_star_loc[localId])
                     );
    }

    barrier(CLK_LOCAL_MEM_FENCE)
    /* Now reduce the sum */

}
