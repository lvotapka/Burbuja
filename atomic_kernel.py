import cupy as cp

kernel_code = r'''
extern "C" __global__
void atomic_add_kernel(const float* masses, const int* indices, float* mass_array, int N) {
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < N) {
        atomicAdd(&mass_array[indices[idx]], masses[idx]);
    }
}
'''

# Compile kernel
module = cp.RawModule(code=kernel_code)
atomic_add_kernel = module.get_function('atomic_add_kernel')

def atomic_add(mass_array, indices, masses):
    N = indices.size
    threads_per_block = 256
    blocks = (N + threads_per_block - 1) // threads_per_block
    atomic_add_kernel((blocks,), (threads_per_block,), 
                      (masses, indices, mass_array, N))