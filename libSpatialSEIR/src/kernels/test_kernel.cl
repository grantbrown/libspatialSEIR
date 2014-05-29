__kernel void test_kernel(__global float* a,
                          __global float* b,
                          __global float* output)
{
    size_t idx = get_global_id(0);
    output[idx] = (a[idx])+(b[idx]);
}
