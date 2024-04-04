/*! @file BCPeriodic_GPU.cu
    @author Youngdae Kim
    @brief GPU implementations of periodic boundary conditions
*/
#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>
#include <mpivars.h>
#include <boundaryconditions.h>

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel of  gpuBCPeriodicU() */
__global__
void BCPeriodicU_kernel(
  int npoints_bounds,
  int npoints_local_wghosts,
  int face,
  int ndims,
  int dim,
  int ghosts,
  int nvars,
  const int * __restrict__ bounds,
  const int * __restrict__ size,
  const int * __restrict__ boundary_is,
  double * __restrict__ phi
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_bounds) {
    int p1 = 0, p2 = 0;
    int index1[GPU_MAX_NDIMS], index2[GPU_MAX_NDIMS];

    _ArrayIndexnD_(ndims,p,bounds,index1,0);
    _ArrayCopy1D_(index1,index2,ndims);
    if (face == 1) {
      index2[dim] = index1[dim] + size[dim]-ghosts;
      _ArrayIndex1DWO_(ndims,size,index1,boundary_is,ghosts,p1);
      _ArrayIndex1D_(ndims,size,index2,ghosts,p2);
    } else if (face == -1) {
      _ArrayIndex1DWO_(ndims,size,index1,boundary_is,ghosts,p1);
      _ArrayIndex1D_(ndims,size,index1,ghosts,p2);
    }
    _ArrayCopy1D_((phi+nvars*p2),(phi+nvars*p1),nvars);
  }
  return;
}

/*! Applies periodic boundary conditions: Implemented by copying the solution
    from the other end of the domain into the physical boundary ghost points.
    \n\n
    **Note**: This function only acts if the the number of processors is 1 along
    the spatial dimension this boundary corresponds to. If there are more than 1
    processors along this dimension, periodicity is handled by MPIExchangeBoundariesnD()
    to minimize communication.
    \sa BCPeriodicU() */
extern "C" int gpuBCPeriodicU(
    void   * __restrict__ b,  /*!< Boundary object of type #DomainBoundary */
    void   * __restrict__ m,  /*!< MPI object of type #MPIVariables */
    int    ndims,             /*!< Number of spatial dimensions */
    int    nvars,             /*!< Number of variables/DoFs per grid point */
    int    * __restrict__ size, /*!<  Integer array with the number of grid points in
                                      each spatial dimensions */
    int    ghosts,              /*!< Number of ghost points */
    double * __restrict__ phi,  /*!< The solution array on which to apply the boundary condition */
    double waqt                 /*!< Current solution time */
)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  MPIVariables   *mpi      = (MPIVariables*)   m;

  int dim   = boundary->dim;
  int face  = boundary->face;

  if ((boundary->on_this_proc) && (mpi->iproc[dim] == 1)) {
    int nblocks;
    nblocks = (boundary->gpu_npoints_bounds-1) / GPU_THREADS_PER_BLOCK + 1;

#if defined(GPU_STAT)
    cudaEvent_t startEvent, stopEvent;
    float milliseconds = 0;

    int memory_accessed = 2*boundary->gpu_npoints_bounds*nvars*sizeof(double);

    checkCuda( cudaEventCreate(&startEvent));
    checkCuda( cudaEventCreate(&stopEvent));

    checkCuda( cudaEventRecord(startEvent, 0) );
#endif

    BCPeriodicU_kernel<<<nblocks,GPU_THREADS_PER_BLOCK>>>(boundary->gpu_npoints_bounds,
      boundary->gpu_npoints_local_wghosts, face, ndims, dim, ghosts, nvars,
      boundary->gpu_bounds, size, boundary->gpu_is, phi
    );
    cudaDeviceSynchronize();

#if defined(GPU_STAT)
    checkCuda( cudaEventRecord(stopEvent, 0) );
    checkCuda( cudaEventSynchronize(stopEvent) );
    checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

    printf("%-50s GPU time (secs) = %.6f bandwidth (GB/s) = %6.2f\n",
           "BCPeriodicU", milliseconds*1e-3,
           (1e-6*(memory_accessed)/milliseconds));
#endif
  }

  return(0);
}

#else

/*! Kernel of  gpuBCPeriodicU() */
__global__
void BCPeriodicU_kernel(
  int npoints_bounds,
  int npoints_local_wghosts,
  int face,
  int ndims,
  int dim,
  int ghosts,
  int nvars,
  const int * __restrict__ bounds,
  const int * __restrict__ size,
  const int * __restrict__ boundary_is,
  double * __restrict__ phi
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_bounds) {
    int p1 = 0, p2 = 0;
    int index1[GPU_MAX_NDIMS], index2[GPU_MAX_NDIMS];

    _ArrayIndexnD_(ndims,p,bounds,index1,0);
    _ArrayCopy1D_(index1,index2,ndims);
    if (face == 1) {
      index2[dim] = index1[dim] + size[dim]-ghosts;
      _ArrayIndex1DWO_(ndims,size,index1,boundary_is,ghosts,p1);
      _ArrayIndex1D_(ndims,size,index2,ghosts,p2);
    } else if (face == -1) {
      _ArrayIndex1DWO_(ndims,size,index1,boundary_is,ghosts,p1);
      _ArrayIndex1D_(ndims,size,index1,ghosts,p2);
    }

    for (int j=0; j<nvars; j++) {
      phi[p1] = phi[p2];
      p1 += npoints_local_wghosts;
      p2 += npoints_local_wghosts;
    }
  }
  return;
}

/*! Applies periodic boundary conditions: Implemented by copying the solution
    from the other end of the domain into the physical boundary ghost points.
    \n\n
    **Note**: This function only acts if the the number of processors is 1 along
    the spatial dimension this boundary corresponds to. If there are more than 1
    processors along this dimension, periodicity is handled by MPIExchangeBoundariesnD()
    to minimize communication.
    \sa BCPeriodicU(), gpuBCPeriodicU() */
extern "C" int gpuBCPeriodicU(
    void   * __restrict__ b,  /*!< Boundary object of type #DomainBoundary */
    void   * __restrict__ m,  /*!< MPI object of type #MPIVariables */
    int    ndims,             /*!< Number of spatial dimensions */
    int    nvars,             /*!< Number of variables/DoFs per grid point */
    int    * __restrict__ size, /*!<  Integer array with the number of grid points in
                                      each spatial dimensions */
    int    ghosts,              /*!< Number of ghost points */
    double * __restrict__ phi,  /*!< The solution array on which to apply the boundary condition */
    double waqt                 /*!< Current solution time */
)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  MPIVariables   *mpi      = (MPIVariables*)   m;

  int dim   = boundary->dim;
  int face  = boundary->face;

  if ((boundary->on_this_proc) && (mpi->iproc[dim] == 1)) {
    int nblocks;
    nblocks = (boundary->gpu_npoints_bounds-1) / GPU_THREADS_PER_BLOCK + 1;

#if defined(GPU_STAT)
    cudaEvent_t startEvent, stopEvent;
    float milliseconds = 0;

    int memory_accessed = 2*boundary->gpu_npoints_bounds*nvars*sizeof(double);

    checkCuda(cudaEventCreate(&startEvent));
    checkCuda(cudaEventCreate(&stopEvent));

    checkCuda(cudaEventRecord(startEvent, 0));
#endif

    BCPeriodicU_kernel<<<nblocks,GPU_THREADS_PER_BLOCK>>>(boundary->gpu_npoints_bounds,
      boundary->gpu_npoints_local_wghosts, face, ndims, dim, ghosts, nvars,
      boundary->gpu_bounds, size, boundary->gpu_is, phi
    );
    cudaDeviceSynchronize();

#if defined(GPU_STAT)
    checkCuda(cudaEventRecord(stopEvent, 0));
    checkCuda(cudaEventSynchronize(stopEvent));
    checkCuda(cudaEventElapsedTime(&milliseconds, startEvent, stopEvent));

    printf("%-50s GPU time (secs) = %.6f bandwidth (GB/s) = %6.2f\n",
           "BCPeriodicU2", milliseconds*1e-3,
           (1e-6*(memory_accessed)/milliseconds));
#endif
  }

  return(0);
}

#endif
