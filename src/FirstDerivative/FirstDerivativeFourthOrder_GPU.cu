/*! @file FirstDerivativeFourthOrder_GPU.cu
    @author Youngdae Kim
    @brief GPU implementation of fourth order finite-difference approximation to the first derivative
*/
#include <basic_gpu.h>
#include <firstderivative.h>
#include <arrayfunctions.h>

#include <mpivars.h>
#include <hypar.h>
typedef MPIVariables  MPIContext;
typedef HyPar         SolverContext;

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel for gpuFirstDerivativeFourthOrderCentral() */
__global__
void FirstDerivativeFourthOrderCentral_boundary_kernel(
    int N_outer,
    int ghosts,
    int ndims,
    int nvars,
    int dir,
    const int *dim,
    const double *f,
    double *Df
)
{
    int j = threadIdx.x + (blockDim.x * blockIdx.x);
    if (j < N_outer) {
        double one_twelve = 1.0/12.0;

        int i, v;
        int indexC[GPU_MAX_NDIMS], index_outer[GPU_MAX_NDIMS], bounds_outer[GPU_MAX_NDIMS];
        _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
        _ArrayIndexnD_(ndims,j,bounds_outer,index_outer,0);
        _ArrayCopy1D_(index_outer,indexC,ndims);

        /* left boundary */
        for (i = -ghosts; i < -ghosts+1; i++) {
            int     qC, qp1, qp2, qp3, qp4;
            indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
            indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
            indexC[dir] = i+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
            indexC[dir] = i+3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp3);
            indexC[dir] = i+4; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp4);
            for (v=0; v<nvars; v++)
              Df[qC*nvars+v] = (-25*f[qC*nvars+v]+48*f[qp1*nvars+v]-36*f[qp2*nvars+v]+16*f[qp3*nvars+v]-3*f[qp4*nvars+v])*one_twelve;
        }
        for (i = -ghosts+1; i < -ghosts+2; i++) {
            int qC, qm1, qp1, qp2, qp3;
            indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
            indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
            indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
            indexC[dir] = i+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
            indexC[dir] = i+3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp3);
            for (v=0; v<nvars; v++)
              Df[qC*nvars+v] = (-3*f[qm1*nvars+v]-10*f[qC*nvars+v]+18*f[qp1*nvars+v]-6*f[qp2*nvars+v]+f[qp3*nvars+v])*one_twelve;
          }
        /* right boundary */
        for (i = dim[dir]+ghosts-2; i < dim[dir]+ghosts-1; i++) {
            int qC, qm3, qm2, qm1, qp1;
            indexC[dir] = i-3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
            indexC[dir] = i-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
            indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
            indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
            indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
            for (v=0; v<nvars; v++)
                Df[qC*nvars+v] = (-f[qm3*nvars+v]+6*f[qm2*nvars+v]-18*f[qm1*nvars+v]+10*f[qC*nvars+v]+3*f[qp1*nvars+v])*one_twelve;
        }
        for (i = dim[dir]+ghosts-1; i < dim[dir]+ghosts; i++) {
            int qC, qm4, qm3, qm2, qm1;
            indexC[dir] = i-4; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm4);
            indexC[dir] = i-3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
            indexC[dir] = i-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
            indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
            indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
            for (v=0; v<nvars; v++)
                Df[qC*nvars+v] = (3*f[qm4*nvars+v]-16*f[qm3*nvars+v]+36*f[qm2*nvars+v]-48*f[qm1*nvars+v]+25*f[qC*nvars+v])*one_twelve;
        }
    }

    return;
}

/*! Kernel for gpuFirstDerivativeFourthOrderCentral() */
__global__
void FirstDerivativeFourthOrderCentral_interior_kernel(
    int ngrid_points,
    int ghosts,
    int ndims,
    int nvars,
    int dir,
    const int *dim,
    const double *f,
    double *Df
)
{
    int i = threadIdx.x + (blockDim.x * blockIdx.x);
    if (i < ngrid_points) {
        /* interior */
        double one_twelve = 1.0/12.0;

        int j, v;
        int qC, qm1, qm2, qp1, qp2;
        int indexC[GPU_MAX_NDIMS], index_outer[GPU_MAX_NDIMS], bounds_outer[GPU_MAX_NDIMS];

        j = i/(dim[dir] + 2*ghosts - 4); /* Outer index */
        _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
        _ArrayIndexnD_(ndims,j,bounds_outer,index_outer,0);
        _ArrayCopy1D_(index_outer,indexC,ndims);

        //i += (-ghosts + 2);
        i = (i % (dim[dir] + 2*ghosts - 4)) + (-ghosts + 2);
        indexC[dir] = i-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
        indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = i+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
        for (v=0; v<nvars; v++)
            Df[qC*nvars+v] = (f[qm2*nvars+v]-8*f[qm1*nvars+v]+8*f[qp1*nvars+v]-f[qp2*nvars+v])*one_twelve;
    }
    return;
}

/*! Computes the fourth-order finite-difference approximation to the first derivative
    (\b Note: not divided by the grid spacing):
    \f{equation}{
      \left(\partial f\right)_i = \left\{ \begin{array}{ll} \frac{1}{12}\left(-25f_i+48f_{i+1}-36f_{i+2}+16f_{i+3}-3f_{i+4}\right) & i=-g \\ \frac{1}{12}\left(-3f_{i-1}-10f_i+18f_{i+1}-6f_{i+2}+f_{i+3}\right) & i = -g+1 \\ \frac{1}{2}\left( f_{i-2}-8f_{i-1}+8f_{i+1}-f_{i+2} \right) & -g+2 \leq i \leq N+g-3 \\ \frac{1}{12}\left( -f_{i-3}+6f_{i-2}-18f_{i-1}+10f_i+3f_{i+1}\right) & i = N+g-2 \\ \frac{1}{12}\left( 3f_{i-4}-16f_{i-3}+36f_{i-2}-48f_{i-1}+25f_i \right) & i = N+g-1 \end{array}\right.
    \f}
    where \f$i\f$ is the grid index along the spatial dimension of the derivative, \f$g\f$ is the number of ghost points, and \f$N\f$ is the number of grid points (not including the ghost points) in the spatial dimension of the derivative.
    \n\n
    Notes:
    + The first derivative is computed at the grid points or the cell centers.
    + The first derivative is computed at the ghost points too. Thus, biased schemes are used
      at and near the boundaries.
    + \b Df and \b f are 1D arrays containing the function and its computed derivatives on a multi-
      dimensional grid. The derivative along the specified dimension \b dir is computed by looping
      through all grid lines along \b dir.

  \sa FirstDerivativeFourthOrderCentral()
*/
int gpuFirstDerivativeFourthOrderCentral(
    double  *Df,  /*!< Array to hold the computed first derivative (with ghost points) */
    double  *f,   /*!< Array containing the grid point function values whose first
                        derivative is to be computed (with ghost points) */
    int     dir,  /*!< The spatial dimension along which the derivative is computed */
    int     bias, /*!< Forward or backward differencing for non-central
                        finite-difference schemes (-1: backward, 1: forward)*/
    void    *s,   /*!< Solver object of type #SolverContext */
    void    *m    /*!< MPI object of type #MPIContext */
)
{
  SolverContext *solver = (SolverContext*) s;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  if ((!Df) || (!f)) {
    fprintf(stderr, "Error in FirstDerivativeFourthOrder(): input arrays not allocated.\n");
    return(1);
  }

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int bounds_outer[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);
  int nblocks = (N_outer - 1) / GPU_THREADS_PER_BLOCK + 1;

#if defined(GPU_STAT)
  cudaEvent_t startEvent, stopEvent;
  float milliseconds;

  checkCuda( cudaEventCreate(&startEvent) );
  checkCuda( cudaEventCreate(&stopEvent) );

  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  FirstDerivativeFourthOrderCentral_boundary_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      N_outer, ghosts, ndims, nvars, dir, solver->gpu_dim_local, f, Df
  );

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  printf("%-50s GPU time (secs) = %.6f\n",
          "FirstDerivativeFourthOrderCentral_boundary", milliseconds*1e-3);
#endif

  int npoints_grid = N_outer*(dim[dir] + 2*ghosts - 4);
  nblocks = (npoints_grid - 1) / GPU_THREADS_PER_BLOCK + 1;

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  FirstDerivativeFourthOrderCentral_interior_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    npoints_grid, ghosts, ndims, nvars, dir, solver->gpu_dim_local, f, Df
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  printf("%-50s GPU time (secs) = %.6f\n",
          "FirstDerivativeFourthOrderCentral_interior", milliseconds*1e-3);
#endif
  return(0);
}

#else

/*! Kernel for gpuFirstDerivativeFourthOrderCentral() */
__global__
void FirstDerivativeFourthOrderCentral_boundary_kernel(
    int N_outer,
    int npoints_local_wghosts,
    int ghosts,
    int ndims,
    int nvars,
    int dir,
    const int *dim,
    const double *f,
    double *Df
)
{
    int j = threadIdx.x + (blockDim.x * blockIdx.x);
    if (j < N_outer) {
        double one_twelve = 1.0/12.0;

        int i, v;
        int indexC[GPU_MAX_NDIMS], index_outer[GPU_MAX_NDIMS], bounds_outer[GPU_MAX_NDIMS];
        _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
        _ArrayIndexnD_(ndims,j,bounds_outer,index_outer,0);
        _ArrayCopy1D_(index_outer,indexC,ndims);

        /* left boundary */
        for (i = -ghosts; i < -ghosts+1; i++) {
            int     qC, qp1, qp2, qp3, qp4;
            indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
            indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
            indexC[dir] = i+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
            indexC[dir] = i+3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp3);
            indexC[dir] = i+4; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp4);
            for (v=0; v<nvars; v++) {
              Df[qC+v*npoints_local_wghosts] = (-25*f[qC+v*npoints_local_wghosts]+48*f[qp1+v*npoints_local_wghosts]-36*f[qp2+v*npoints_local_wghosts]+16*f[qp3+v*npoints_local_wghosts]-3*f[qp4+v*npoints_local_wghosts])*one_twelve;
            }
        }
        for (i = -ghosts+1; i < -ghosts+2; i++) {
            int qC, qm1, qp1, qp2, qp3;
            indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
            indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
            indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
            indexC[dir] = i+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
            indexC[dir] = i+3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp3);
            for (v=0; v<nvars; v++)
              Df[qC+v*npoints_local_wghosts] = (-3*f[qm1+v*npoints_local_wghosts]-10*f[qC+v*npoints_local_wghosts]+18*f[qp1+v*npoints_local_wghosts]-6*f[qp2+v*npoints_local_wghosts]+f[qp3+v*npoints_local_wghosts])*one_twelve;
          }
        /* right boundary */
        for (i = dim[dir]+ghosts-2; i < dim[dir]+ghosts-1; i++) {
            int qC, qm3, qm2, qm1, qp1;
            indexC[dir] = i-3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
            indexC[dir] = i-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
            indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
            indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
            indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
            for (v=0; v<nvars; v++)
                Df[qC+v*npoints_local_wghosts] = (-f[qm3+v*npoints_local_wghosts]+6*f[qm2+v*npoints_local_wghosts]-18*f[qm1+v*npoints_local_wghosts]+10*f[qC+v*npoints_local_wghosts]+3*f[qp1+v*npoints_local_wghosts])*one_twelve;
        }
        for (i = dim[dir]+ghosts-1; i < dim[dir]+ghosts; i++) {
            int qC, qm4, qm3, qm2, qm1;
            indexC[dir] = i-4; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm4);
            indexC[dir] = i-3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
            indexC[dir] = i-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
            indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
            indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
            for (v=0; v<nvars; v++)
                Df[qC+v*npoints_local_wghosts] = (3*f[qm4+v*npoints_local_wghosts]-16*f[qm3+v*npoints_local_wghosts]+36*f[qm2+v*npoints_local_wghosts]-48*f[qm1+v*npoints_local_wghosts]+25*f[qC+v*npoints_local_wghosts])*one_twelve;
        }
    }

    return;
}

/*! Kernel for gpuFirstDerivativeFourthOrderCentral() */
__global__
void FirstDerivativeFourthOrderCentral_interior_kernel(
    int ngrid_points,
    int npoints_local_wghosts,
    int ghosts,
    int ndims,
    int nvars,
    int dir,
    const int *dim,
    const double *f,
    double *Df
)
{
    int i = threadIdx.x + (blockDim.x * blockIdx.x);
    if (i < ngrid_points) {
        /* interior */
        double one_twelve = 1.0/12.0;

        int j, v;
        int qC, qm1, qm2, qp1, qp2;
        int indexC[GPU_MAX_NDIMS], index_outer[GPU_MAX_NDIMS], bounds_outer[GPU_MAX_NDIMS];

        j = i/(dim[dir] + 2*ghosts - 4); /* Outer index */
        _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
        _ArrayIndexnD_(ndims,j,bounds_outer,index_outer,0);
        _ArrayCopy1D_(index_outer,indexC,ndims);

        //i += (-ghosts + 2);
        i = (i % (dim[dir] + 2*ghosts - 4)) + (-ghosts + 2);
        indexC[dir] = i-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
        indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = i+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
        for (v=0; v<nvars; v++)
            Df[qC+v*npoints_local_wghosts] = (f[qm2+v*npoints_local_wghosts]-8*f[qm1+v*npoints_local_wghosts]+8*f[qp1+v*npoints_local_wghosts]-f[qp2+v*npoints_local_wghosts])*one_twelve;
    }
    return;
}

/*! Computes the fourth-order finite-difference approximation to the first derivative
    (\b Note: not divided by the grid spacing):
    \f{equation}{
      \left(\partial f\right)_i = \left\{ \begin{array}{ll} \frac{1}{12}\left(-25f_i+48f_{i+1}-36f_{i+2}+16f_{i+3}-3f_{i+4}\right) & i=-g \\ \frac{1}{12}\left(-3f_{i-1}-10f_i+18f_{i+1}-6f_{i+2}+f_{i+3}\right) & i = -g+1 \\ \frac{1}{2}\left( f_{i-2}-8f_{i-1}+8f_{i+1}-f_{i+2} \right) & -g+2 \leq i \leq N+g-3 \\ \frac{1}{12}\left( -f_{i-3}+6f_{i-2}-18f_{i-1}+10f_i+3f_{i+1}\right) & i = N+g-2 \\ \frac{1}{12}\left( 3f_{i-4}-16f_{i-3}+36f_{i-2}-48f_{i-1}+25f_i \right) & i = N+g-1 \end{array}\right.
    \f}
    where \f$i\f$ is the grid index along the spatial dimension of the derivative, \f$g\f$ is the number of ghost points, and \f$N\f$ is the number of grid points (not including the ghost points) in the spatial dimension of the derivative.
    \n\n
    Notes:
    + The first derivative is computed at the grid points or the cell centers.
    + The first derivative is computed at the ghost points too. Thus, biased schemes are used
      at and near the boundaries.
    + \b Df and \b f are 1D arrays containing the function and its computed derivatives on a multi-
      dimensional grid. The derivative along the specified dimension \b dir is computed by looping
      through all grid lines along \b dir.

  \sa FirstDerivativeFourthOrderCentral(), gpuFirstDerivativeFourthOrderCentral()
*/
int gpuFirstDerivativeFourthOrderCentral(
    double  *Df,  /*!< Array to hold the computed first derivative (with ghost points) */
    double  *f,   /*!< Array containing the grid point function values whose first
                        derivative is to be computed (with ghost points) */
    int     dir,  /*!< The spatial dimension along which the derivative is computed */
    int     bias, /*!< Forward or backward differencing for non-central
                        finite-difference schemes (-1: backward, 1: forward)*/
    void    *s,   /*!< Solver object of type #SolverContext */
    void    *m    /*!< MPI object of type #MPIContext */
)
{
  SolverContext *solver = (SolverContext*) s;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  if ((!Df) || (!f)) {
    fprintf(stderr, "Error in FirstDerivativeFourthOrder(): input arrays not allocated.\n");
    return(1);
  }

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int bounds_outer[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);
  int nblocks = (N_outer - 1) / GPU_THREADS_PER_BLOCK + 1;

#if defined(GPU_STAT)
  cudaEvent_t startEvent, stopEvent;
  float milliseconds;

  checkCuda( cudaEventCreate(&startEvent) );
  checkCuda( cudaEventCreate(&stopEvent) );

  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  FirstDerivativeFourthOrderCentral_boundary_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      N_outer, solver->npoints_local_wghosts, ghosts, ndims, nvars, dir, solver->gpu_dim_local, f, Df
  );

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  printf("%-50s GPU time (secs) = %.6f\n",
          "FirstDerivativeFourthOrderCentral_boundary", milliseconds*1e-3);
#endif

  int npoints_grid = N_outer*(dim[dir] + 2*ghosts - 4);
  nblocks = (npoints_grid - 1) / GPU_THREADS_PER_BLOCK + 1;

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  FirstDerivativeFourthOrderCentral_interior_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    npoints_grid, solver->npoints_local_wghosts, ghosts, ndims, nvars, dir, solver->gpu_dim_local, f, Df
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  printf("%-50s GPU time (secs) = %.6f\n",
          "FirstDerivativeFourthOrderCentral_interior", milliseconds*1e-3);
#endif

  return(0);
}

#endif
