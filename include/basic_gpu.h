#ifndef _BASIC_GPU_H_
#define _BASIC_GPU_H_

#include <basic.h>

#define GPU_ONE_SIXTH           (1.0/6.0)
#define GPU_THREADS_PER_BLOCK   64
#define GPU_MAX_NDIMS            3
#define GPU_MAX_NVARS            5

#if defined(HAVE_CUDA) && defined(__NVCC__)
    inline
    cudaError_t checkCuda(cudaError_t result)
    {
    #if defined(DEBUG) || defined(_DEBUG)
        if (result != cudaSuccess) {
            fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
            assert(result == cudaSuccess);
        }
    #endif
        return result;
    }
#endif

#endif /* _BASIC_GPU_H_ */
