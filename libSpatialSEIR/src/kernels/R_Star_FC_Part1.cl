#pragma OPENCL EXTENSION cl_khr_fp64 : enable
__kernel void FC_R_Star_Part1(__global int* R_star,
                              __global int* S_star,
                              __global int* R,
                              __global int* I,
                              __global double p_rs,
                              __global double p_ir,
                              __global double output)
{
    size_t i = get_global_id(0);


    double p_ir_val = *p_ir;
    double ln_p_ir = log(p_ir_val);
    double ln_1m_p_ir = log(1-p_ir_val);
    int Rstar_val = R_star[i]; 
    int Sstar_val = S_star[i]; 
    int I_val = I[i]; 
    int R_val = R[i]; 
    double p_rs_val = p_rs[i];

    if (Rstar_val < 0 || Rstar_val > I_val || 
                  Sstar_val > R_val)
    {
        output[i] = -INFINITY;
        return;
    }
    
    output[i] = (ln_p_ir*Rstar_val + 
                 ln_1m_p_ir*(I_val - Rstar_val) + 
                 log(p_rs_val)*Sstar_val + log(1-p_rs_val)*(R_val) + 
                 R_val*log(1.0*R_val) - 
                 (R_val-Sstar_val)*log(1.0*R_val-Sstar_val) - Sstar_val +  
                 I_val*log(1.0*I_val) - Rstar_val*log(1.0*Rstar_val) - 
                 (I_val - Rstar_val)*log(1.0*I_val - Rstar_val)
                 );
}
