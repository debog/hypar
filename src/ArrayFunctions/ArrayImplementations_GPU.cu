/*! @file ArrayImplementations_GPU.cu
    @author Youngdae Kim
    @brief GPU implementations of array functions.
*/

#include <assert.h>
#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>

/*! Element-wise copy \a y = \a x, where \a x, \a y are 1-dimensional arrays
    of length \a size. \sa #_ArrayCopy1D_ */
__global__
void ArrayCopy1D_kernel(  const double* x,  /*!< copy-from array*/
                          double*       y,  /*!< copy-to array */
                          int n             /*!< size of array */ )
{
    int tx = threadIdx.x + (blockIdx.x * blockDim.x);
    if (tx < n) y[tx] = x[tx];
    return;
}

/*! Set all elements of a 1-dimensional array \a x (any datatype)
    of length \a size to a scalar \a value. \sa #_ArraySetValue_ */
__global__
void ArraySetValue_kernel(  double* x,      /*!< array*/
                            int     n,      /*!< size of array */
                            double  value   /*!< scalar value */ )
{
    int tx = threadIdx.x + (blockIdx.x * blockDim.x);
    if (tx < n) x[tx] = value;
    return;
}

/*! \sa #_ArrayAXPY_

  Element-wise AXPY \a y = \a a \a x + \a y,
  where \a a is a scalar, and \a x, \a y, \a z are
  1-dimensional arrays of length \a size. */
__global__
void ArrayAXPY_kernel(  const double* x,  /*!< x */
                        double        a,  /*!< a */
                        double*       y,  /*!< y */
                        int           n   /*!< size of array */ )
{
    int tx = threadIdx.x + (blockIdx.x * blockDim.x);
    if (tx < n) y[tx] += a*x[tx];
    return;
}

/*! \sa #_ArrayBlockMultiply_

  Given two arrays: \a x of size \a n*bs, and \a a of size \a n, this function
  implements: \a x[i][j] *= \a a[i] where \a i = 1,...,\a n, j = 1,...,\a bs,
  and \a x is stored as a 1D array in row-major format, i.e., \a x[i][j] = \a x[i*bs+j]. */
__global__
void ArrayBlockMultiply_kernel( double*       x,  /*!< x */
                                const double* a,  /*!< a */
                                int           n,  /*!< size of array */
                                int           bs  /*!< block size */)
{
    int tx = threadIdx.x + (blockIdx.x * blockDim.x);
    if (tx < n) {
        for (int i = 0; i < bs; i++) x[tx*bs + i] *= a[tx];
    }
}

/*! Alternative implementation of #_ArrayCopy1D_ */
__global__
void ArrayCopy1DNewScheme_kernel( const double* __restrict__ src,   /*!< source array */
                                  double*       __restrict__ dest,  /*!< destination array */
                                  int           npoints,            /*!< number of points */
                                  int           nvars               /*!< number of components */ )
{
    int p = blockDim.x * blockIdx.x + threadIdx.x;
    if (p < npoints) {
        for (int v=0; v<nvars; v++) {
            dest[p+v*npoints] = src[p*nvars+v];
        }
    }
    return;
}

/*! Set device */
void gpuSetDevice(int device /*!< device */)
{
    cudaSetDevice(device);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf(stderr,"Error in gpuSetDevice(): device=%d error message=\"%s\"\n", device, cudaGetErrorString(err));
    }
}

/*! GPU memory copy */
void gpuMemcpy( void*               dest, /*!< destination */
                const void*         src,  /*!< source */
                size_t              count,/*!< count */
                enum gpuMemcpyKind  kind  /*!< kind of copy */ )
{
    switch (kind) {
        case gpuMemcpyHostToDevice:
            checkCuda( cudaMemcpy(dest, src, count, cudaMemcpyHostToDevice) );
            break;
        case gpuMemcpyDeviceToHost:
            checkCuda( cudaMemcpy(dest, src, count, cudaMemcpyDeviceToHost) );
            break;
        case gpuMemcpyDeviceToDevice:
            checkCuda( cudaMemcpy(dest, src, count, cudaMemcpyDeviceToDevice) );
            break;
        default:
            fprintf(stderr, "Error: invalid gpuMemcpyKind: %d\n", kind);
            assert(0);
            break;
    }
    return;
}

/*! Allocate memory */
void gpuMalloc( void**  devPtr, /*!< pointer to memory */
                size_t  size    /*!< size of memory */ )
{
    cudaMalloc(devPtr, size);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf(  stderr,"Error in gpuMalloc(): size=%d, error message=\"%s\"\n", size,
                cudaGetErrorString(err) );
    }
    return;
}

/*! Set value */
void gpuMemset( void*   devPtr, /*!< Pointer to memory */
                int     value,  /*!< value to set */
                size_t  count   /*!< size of data */ )
{
    checkCuda( cudaMemset(devPtr, value, count) );
    return;
}

/*! deallocate memory */
void gpuFree(void *devPtr /*!< Pointer to memory */)
{
    checkCuda( cudaFree(devPtr) );
    return;
}

/*! Element-wise copy \a y = \a x, where \a x, \a y are 1-dimensional arrays
    of length \a size. \sa #_ArrayCopy1D_, ArrayCopy1D_kernel() */
void gpuArrayCopy1D(  const double* x,  /*!< copy-from array*/
                      double*       y,  /*!< copy-to array */
                      int           n   /*!< size of array */ )
{
    int nblocks = (n - 1) / GPU_THREADS_PER_BLOCK + 1;
    ArrayCopy1D_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(x, y, n);
    cudaDeviceSynchronize();
    return;
}

/*! Set all elements of a 1-dimensional array \a x (any datatype)
    of length \a size to a scalar \a value.
    \sa #_ArraySetValue_, ArraySetValue_kernel() */
void gpuArraySetValue(  double* devPtr, /*!< array*/
                        int     n,      /*!< size of array */
                        double  value   /*!< scalar value */ )
{
    int nblocks = (n - 1) / GPU_THREADS_PER_BLOCK + 1;
    ArraySetValue_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(devPtr, n, value);
    cudaDeviceSynchronize();
    return;
}

/*! \sa #_ArrayAXPY_, ArrayAXPY_kernel()

  Element-wise AXPY \a y = \a a \a x + \a y,
  where \a a is a scalar, and \a x, \a y, \a z are
  1-dimensional arrays of length \a size.
*/
void gpuArrayAXPY(  const double*   x, /*!< x */
                    double          a, /*!< a */
                    double*         y, /*!< y */
                    int             n  /*!< size of array */ )
{
    int nblocks = (n - 1) / GPU_THREADS_PER_BLOCK + 1;
    ArrayAXPY_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(x, a, y, n);
    cudaDeviceSynchronize();
    return;
}

/*! \sa #_ArrayBlockMultiply_, ArrayBlockMultiply_kernel()

  Given two arrays: \a x of size \a n*bs, and \a a of size \a n, this function
  implements: \a x[i][j] *= \a a[i] where \a i = 1,...,\a n, j = 1,...,\a bs,
  and \a x is stored as a 1D array in row-major format, i.e., \a x[i][j] = \a x[i*bs+j]. */
void gpuArrayBlockMultiply( double*       x,  /*!< x */
                            const double* a,  /*!< a */
                            int           n,  /*!< size of array */
                            int           bs  /*!< block size */ )
{
    int nblocks = (n - 1) / GPU_THREADS_PER_BLOCK + 1;
    ArrayBlockMultiply_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(x, a, n, bs);
    cudaDeviceSynchronize();
    return;
}

/*! Returns the sum-of-squares of the elements in an n-D array (useful for L_2 norm)
    \sa ArraySumSquarenD() */
double gpuArraySumSquarenD(int nvars,  /*!< number of elements at one array location,
                                            can be > 1 for systems of equations */
                           int ndims,  /*!< number of dimensions */
                           int *dim,   /*!< integer array of size in each dimension */
                           int ghosts, /*!< number of ghost points in the array x */
                           int *index, /*!< pre-allocated (by the calling function)
                                            integer array of size ndims */
                           double *x   /*!< the array */
                          )
{
    double sum = 0;
    printf("gpuArraySumSquarenD hasn't been implemented, yet.\n");
    exit(0);
    return (sum);
}

/*! Alternative implementation of #_ArrayCopy1D_ */
void gpuArrayCopy1DNewScheme( const double* src,      /*!< source array */
                              double*       dest,     /*!< destination array */
                              int           npoints,  /*!< number of points */
                              int           nvars     /*!< number of components */ )
{
    int nblocks = (npoints-1) / GPU_THREADS_PER_BLOCK + 1;
    ArrayCopy1DNewScheme_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(src, dest, npoints, nvars);
    cudaDeviceSynchronize();
    return;
}

/*! Check if two arrays are equal, if not, report the difference  */
/*! Check if two arrays are equal, if not, report the difference */
void gpuArrayCheckEqual( const char*   msg,      /*!< message */
                         const double* var_adj,  /*!< array */
                         const double* var_sep,  /*!< array */
                         int           npoints   /*!< size of array */ )
{
    double *h_var_adj = (double *) malloc(npoints*sizeof(double));
    double *h_var_sep = (double *) malloc(npoints*sizeof(double));

    gpuMemcpy(h_var_adj, var_adj, npoints*sizeof(double), gpuMemcpyDeviceToHost);
    gpuMemcpy(h_var_sep, var_sep, npoints*sizeof(double), gpuMemcpyDeviceToHost);

    double max_err = 0.0;
    for (int j=0; j<npoints; j++) {
        if (h_var_sep[j] != h_var_adj[j]) {
            max_err = max(max_err, fabs(h_var_sep[j]-h_var_adj[j]));
        }
    }

    free(h_var_adj);
    free(h_var_sep);

    if (max_err > 1e-10) {
        printf("gpuArrayCheckEqual: %-30s max_err = %e\n", msg, max_err);
        exit(0);
    }
    return;
}

