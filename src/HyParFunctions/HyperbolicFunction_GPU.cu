/*! @file HyperbolicFunction_GPU.cu
    @author Youngdae Kim
    @brief Compute the hyperbolic term of the governing equations
*/

#include <arrayfunctions_gpu.h>
#include <basic_gpu.h>
#include <mpivars.h>
#include <hypar.h>

static int ReconstructHyperbolic (double*,double*,double*,double*,int,void*,void*,double,int,
                                  int(*)(double*,double*,double*,double*,
                                         double*,double*,int,void*,double));
static int DefaultUpwinding      (double*,double*,double*,double*,
                                  double*,double*,int,void*,double);

/*! Kernel for stage boundary integral calculation */
template <int block_size>
__global__
void StageBoundaryIntegral_kernel(
  int n,
  int nvars,
  int sign,
  const double * __restrict__ FluxI,
  double       * __restrict__ StageBoundaryIntegral
)
{
  extern __shared__ double smem[];

  int tid        = threadIdx.x;
  int grid_size  = block_size * gridDim.x;
  int i          = blockIdx.x * block_size + threadIdx.x;
  int j;
  double tid_sum[GPU_MAX_NVARS] = { 0 };

  while (i < n) {
    for (j=0; j<nvars; j++) tid_sum[j] += FluxI[i*nvars+j];
    i += grid_size;
  }
  for (j=0; j<nvars; j++) smem[tid*nvars+j] = tid_sum[j];
  __syncthreads();


  if (block_size >= 512 && tid < 256) {
    for (j=0; j<nvars; j++) smem[tid*nvars+j] = tid_sum[j] = tid_sum[j] + smem[(tid+256)*nvars+j];
  }
  __syncthreads();
  if (block_size >= 256 && tid < 128) {
    for (j=0; j<nvars; j++) smem[tid*nvars+j] = tid_sum[j] = tid_sum[j] + smem[(tid+128)*nvars+j];
  }
  __syncthreads();
  if (block_size >= 128 && tid < 64) {
    for (j=0; j<nvars; j++) smem[tid*nvars+j] = tid_sum[j] = tid_sum[j] + smem[(tid+64)*nvars+j];
  }
  __syncthreads();

  if (tid < 32) {
    if (block_size >= 64) {
      for (j=0; j<nvars; j++) tid_sum[j] += smem[(tid+32)*nvars+j];
    }
    for (int offset=16; offset>0; offset/=2) {
      for (j=0; j<nvars; j++) {
        tid_sum[j] += __shfl_down_sync(0xffffffff, tid_sum[j], offset);
      }
    }
  }
  __syncthreads();

  if (tid == 0) {
    for (j=0; j<nvars; j++) StageBoundaryIntegral[j] = (sign == 1) ? tid_sum[j] : -tid_sum[j];
  }
  return;
}

#ifdef CUDA_VAR_ORDERDING_AOS

/*! 3D Kernel for gpuHyperbolicFunction() */
__global__
void HyperbolicFunction_dim3_kernel(
    int npoints_grid,
    int d,
    int nvars,
    int ghosts,
    int offset,
    const int    * __restrict__ dim,
    const double * __restrict__ FluxI,
    const double * __restrict__ dxinv,
    double       * __restrict__ hyp,
    double       * __restrict__ FluxI_p1,
    double       * __restrict__ FluxI_p2
)
{
    int tx = threadIdx.x + (blockDim.x * blockIdx.x);
    if (tx < npoints_grid) {
        int v, p, p1, p2, b;
        int dim_interface[3];
        int index[3], index1[3], index2[3];

        _ArrayIndexnD_(3,tx,dim,index,0);
        _ArrayCopy1D_(dim,dim_interface,3);
        dim_interface[d] += 1;

        _ArrayCopy1D_(index,index1,3);
        _ArrayCopy1D_(index,index2,3); index2[d]++;
        _ArrayIndex1D_(3,dim          ,index ,ghosts,p);
        _ArrayIndex1D_(3,dim_interface,index1,0     ,p1);
        _ArrayIndex1D_(3,dim_interface,index2,0     ,p2);
        for (v=0; v<nvars; v++) hyp[nvars*p+v] += dxinv[offset+ghosts+index[d]]
                                                * (FluxI[nvars*p2+v]-FluxI[nvars*p1+v]);
        if (index[d] == 0 || index[d] == dim[d]-1) {
          if (d == 0) {
            b = index[1] + index[2]*dim[1];
          } else if (d == 1) {
            b = index[0] + index[2]*dim[0];
          } else {
            b = index[0] + index[1]*dim[0];
          }
          if (index[d] == 0)
            for (v=0; v<nvars; v++) FluxI_p1[b*nvars+v] = FluxI[nvars*p1+v];
          else
            for (v=0; v<nvars; v++) FluxI_p2[b*nvars+v] = FluxI[nvars*p2+v];
        }
    }
    return;
}

/*! 2D Kernel for gpuHyperbolicFunction() */
__global__
void HyperbolicFunction_dim2_kernel(
    int npoints_grid,
    int d,
    int nvars,
    int ghosts,
    int offset,
    const int    * __restrict__ dim,
    const double * __restrict__ FluxI,
    const double * __restrict__ dxinv,
    double       * __restrict__ hyp,
    double       * __restrict__ FluxI_p1,
    double       * __restrict__ FluxI_p2
)
{
    int tx = threadIdx.x + (blockDim.x * blockIdx.x);
    if (tx < npoints_grid) {
        int v, p, p1, p2, b;
        int dim_interface[2];
        int index[2], index1[2], index2[2];

        _ArrayIndexnD_(2,tx,dim,index,0);
        _ArrayCopy1D_(dim,dim_interface,2);
        dim_interface[d] += 1;

        _ArrayCopy1D_(index,index1,2);
        _ArrayCopy1D_(index,index2,2); index2[d]++;
        _ArrayIndex1D_(2,dim          ,index ,ghosts,p);
        _ArrayIndex1D_(2,dim_interface,index1,0     ,p1);
        _ArrayIndex1D_(2,dim_interface,index2,0     ,p2);
        for (v=0; v<nvars; v++) hyp[nvars*p+v] += dxinv[offset+ghosts+index[d]]
                                                * (FluxI[nvars*p2+v]-FluxI[nvars*p1+v]);
        if (index[d] == 0 || index[d] == dim[d]-1) {
          b = (d == 0) ? index[1] : index[0];
          if (index[d] == 0)
            for (v=0; v<nvars; v++) FluxI_p1[b*nvars+v] = FluxI[nvars*p1+v];
          else
            for (v=0; v<nvars; v++) FluxI_p2[b*nvars+v] = FluxI[nvars*p2+v];
        }
    }
    return;
}

/*! This function computes the hyperbolic term of the governing equations, expressed as follows:-
    \f{equation}{
      {\bf F} \left({\bf u}\right) = \sum_{d=0}^{D-1} \frac {\partial {\bf f}_d\left({\bf u}\right)} {\partial x_d}
    \f}
    using a conservative finite-difference discretization is space:
    \f{equation}{
      \hat{\bf F} \left({\bf u}\right) = \sum_{d=0}^{D-1} \frac {1} {\Delta x_d} \left[ \hat{\bf f}_{d,j+1/2} - \hat{\bf f}_{d,j-1/2} \right]
    \f}
    where \f$d\f$ denotes the spatial dimension, \f$D\f$ denotes the total number of spatial dimensions, the hat denotes
    the discretized quantity, and \f$j\f$ is the grid coordinate along dimension \f$d\f$.
    The approximation to the flux function \f${\bf f}_d\f$ at the interfaces \f$j\pm1/2\f$, denoted by \f$\hat{\bf f}_{d,j\pm 1/2}\f$,
    are computed using the function ReconstructHyperbolic().
*/
extern "C" int gpuHyperbolicFunction(
    double  *hyp, /*!< Array to hold the computed hyperbolic term (shares the same layout as u */
    double  *gpu_u,   /*!< Solution array */
    void    *s,      /*!< Solver object of type #HyPar */
    void    *m,      /*!< MPI object of type #MPIVariables */
    double  t,       /*!< Current simulation time */
    int     LimFlag, /*!< Flag to indicate if the nonlinear coefficients for solution-dependent
                            interpolation method should be recomputed (see ReconstructHyperbolic() for
                            an explanation on why this is needed) */
    /*! Function pointer to the flux function for the hyperbolic term */
    int(*FluxFunction)(double*,double*,int,void*,double),
    /*! Function pointer to the upwinding function for the hyperbolic term */
    int(*UpwindFunction)(double*,double*,double*,double*,double*,double*,int,void*,double)
)
{
    HyPar         *solver = (HyPar*)        s;
    MPIVariables  *mpi    = (MPIVariables*) m;
    int           d;
    double        *gpu_FluxI  = solver->fluxI; /* interface flux     */
    double        *gpu_FluxC  = solver->fluxC; /* cell centered flux */

    int     ndims     = solver->ndims;
    int     nvars     = solver->nvars;
    int     ghosts    = solver->ghosts;
    int     *dim      = solver->dim_local;
    int     size      = solver->npoints_local_wghosts;
    double  *gpu_x     = solver->gpu_x;
    double  *gpu_dxinv = solver->gpu_dxinv;
    double  *gpu_FluxI_p1, *gpu_FluxI_p2;

    LimFlag = (LimFlag && solver->flag_nonlinearinterp && solver->SetInterpLimiterVar);

    gpuMemset(hyp, 0, size*nvars*sizeof(double));
    gpuMemset(solver->StageBoundaryBuffer, 0, solver->StageBoundaryBuffer_size*sizeof(double));
    gpuMemset(solver->StageBoundaryIntegral, 0, 2*ndims*nvars*sizeof(double));

    if (!FluxFunction) return(0); /* zero hyperbolic term */
    /*solver->count_hyp++;*/

    int npoints_grid = 1; for (d = 0; d < ndims; d++) npoints_grid *= dim[d];
    int nblocks = (npoints_grid - 1) / GPU_THREADS_PER_BLOCK + 1;
    int offset = 0;

#if defined(GPU_STAT)
    cudaEvent_t start, stop;
    float milliseconds = 0;

    checkCuda( cudaEventCreate(&start) );
    checkCuda( cudaEventCreate(&stop) );
#endif

    for (d = 0; d < ndims; d++) {
        /* evaluate cell-centered flux */
        FluxFunction(gpu_FluxC,gpu_u,d,solver,t);
        /* compute interface fluxes */
        ReconstructHyperbolic(gpu_FluxI,gpu_FluxC,gpu_u,gpu_x+offset,d,solver,mpi,t,LimFlag,UpwindFunction);

        gpu_FluxI_p1 = solver->StageBoundaryBuffer+(solver->gpu_npoints_boundary_offset[d]*nvars);
        gpu_FluxI_p2 = solver->StageBoundaryBuffer+(solver->gpu_npoints_boundary_offset[d]+solver->gpu_npoints_boundary[d])*nvars;

#if defined(GPU_STAT)
        checkCuda( cudaEventRecord(start, 0) );
#endif

        if (ndims == 3) {
          HyperbolicFunction_dim3_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
              npoints_grid, d, nvars, ghosts, offset,
              solver->gpu_dim_local, gpu_FluxI, gpu_dxinv,
              hyp, gpu_FluxI_p1, gpu_FluxI_p2
          );
        } else if (ndims == 2) {
          HyperbolicFunction_dim2_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
              npoints_grid, d, nvars, ghosts, offset,
              solver->gpu_dim_local, gpu_FluxI, gpu_dxinv,
              hyp, gpu_FluxI_p1, gpu_FluxI_p2
          );
        } else {
          fprintf(stderr,"gpuHyperbolicFunction(): Not implemented for ndims = %d!\n", ndims);
          exit(1);
        }
        cudaDeviceSynchronize();

#if defined(GPU_STAT)
        checkCuda( cudaEventRecord(stop, 0) );
        checkCuda( cudaEventSynchronize(stop) );
        checkCuda( cudaEventElapsedTime(&milliseconds, start, stop) );
        printf("%-50s GPU time = %.6f dir = %d bandwidth (GB/s) = %.2f\n",
                "HyperbolicFunction", milliseconds*1e-3, d, (1e-6*2*npoints_grid*nvars*sizeof(double))/milliseconds);

        checkCuda( cudaEventRecord(start, 0) );
#endif

        StageBoundaryIntegral_kernel<GPU_THREADS_PER_BLOCK><<<1, GPU_THREADS_PER_BLOCK, GPU_THREADS_PER_BLOCK*nvars*sizeof(double)>>>(
          solver->gpu_npoints_boundary[d], nvars, -1, gpu_FluxI_p1,
          solver->StageBoundaryIntegral + 2*d*nvars
        );
        StageBoundaryIntegral_kernel<GPU_THREADS_PER_BLOCK><<<1, GPU_THREADS_PER_BLOCK, GPU_THREADS_PER_BLOCK*nvars*sizeof(double)>>>(
          solver->gpu_npoints_boundary[d], nvars, 1, gpu_FluxI_p2,
          solver->StageBoundaryIntegral + (2*d+1)*nvars
        );
        cudaDeviceSynchronize();

#if defined(GPU_STAT)
        checkCuda( cudaEventRecord(stop, 0) );
        checkCuda( cudaEventSynchronize(stop) );
        checkCuda( cudaEventElapsedTime(&milliseconds, start, stop) );
        printf("%-50s GPU time = %.6f dir = %d bandwidth (GB/s) = %.2f\n",
                "StageBoundaryIntegral", milliseconds*1e-3, d, (1e-6*2*solver->gpu_npoints_boundary[d]*nvars*sizeof(double))/milliseconds);
#endif

        offset += dim[d] + 2*ghosts;
    }

    if (solver->flag_ib) gpuArrayBlockMultiply(hyp, solver->gpu_iblank, size, nvars);

#if defined(GPU_STAT)
    checkCuda(cudaEventDestroy(start));
    checkCuda(cudaEventDestroy(stop));
#endif

    return(0);
}

/*! This function computes the numerical flux \f$\hat{\bf f}_{j+1/2}\f$ at the interface from the cell-centered
    flux function \f${\bf f}_j\f$. This happens in two steps:-

    \b Interpolation: High-order accurate approximations to the flux at the interface \f$j+1/2\f$ are computed from
    the cell-centered flux with left- and right-biased interpolation methods. This is done by the
    #HyPar::InterpolateInterfacesHyp function. This can be expressed as follows:
    \f{align}{
      \hat{\bf f}^L_{j+1/2} &= \mathcal{I}\left({\bf f}_j,+1\right), \\
      \hat{\bf f}^R_{j+1/2} &= \mathcal{I}\left({\bf f}_j,-1\right),
    \f}
    where the \f$\pm 1\f$ indicates the interpolation bias, and \f$\mathcal{I}\f$ is the interpolation operator
    pointed to by #HyPar::InterpolateInterfacesHyp (see \b src/InterpolationFunctions for all the available operators).
    The interface values of the solution are similarly computed:
    \f{align}{
      \hat{\bf u}^L_{j+1/2} &= \mathcal{I}\left({\bf u}_j,+1\right), \\
      \hat{\bf u}^R_{j+1/2} &= \mathcal{I}\left({\bf u}_j,-1\right),
    \f}
    The specific choice of \f$\mathcal{I}\f$ is set based on #HyPar::spatial_scheme_hyp.

    \b Upwinding: The final flux at the interface is computed as
    \f{equation}{
      \hat{\bf f}_{j+1/2} = \mathcal{U}\left( \hat{\bf f}^L_{j+1/2}, \hat{\bf f}^R_{j+1/2}, \hat{\bf u}^L_{j+1/2}, \hat{\bf u}^R_{j+1/2} \right),
    \f}
    where \f$\mathcal{U}\f$ denotes the upwinding function UpwindFunction() passed as an argument (if NULL, DefaultUpwinding() is used). The
    upwinding function is specified by the physical model.

    \b Note:
    Solution-dependent, nonlinear interpolation methods (such as WENO, CRWENO) are implemented in a way that separates the calculation of the
    nonlinear interpolation weights (based on, say, the smoothness of the flux function), and the actual evaluation of the interpolant, into
    different functions. This allows the flexibility to choose if and when the nonlinear coefficients are evaluated (or previously computed
    values are reused). Some possible scenarios are:
    + For explicit time integration, they are computed every time the hyperbolic flux term is being computed.
    + For implicit time integration, consistency or linearization may require that they be computed and "frozen"
      at the beginning of a time step or stage.

    The argument \b LimFlag controls this behavior:
    + LimFlag = 1 means recompute the nonlinear coefficients.
    + LimFlag = 0 means reuse the the previously computed coefficients.
*/
int ReconstructHyperbolic(
  double  *gpu_fluxI,     /*!< Array to hold the computed interface fluxes. This array does not
                           have ghost points. The dimensions are the same as those of u without
                           ghost points in all dimensions, except along dir, where it is one more */
  double  *gpu_fluxC,     /*!< Array of the flux function computed at the cell centers
                           (same layout as u) */
  double  *gpu_u,         /*!< Solution array */
  double  *gpu_x,         /*!< Array of spatial coordinates */
  int     dir,        /*!< Spatial dimension along which to reconstruct the interface fluxes */
  void    *s,         /*!< Solver object of type #HyPar */
  void    *m,         /*!< MPI object of type #MPIVariables */
  double  t,          /*!< Current solution time */
  int     LimFlag,    /*!< Flag to indicate if the nonlinear coefficients for solution-dependent
                           interpolation method should be recomputed */
  /*! Function pointer to the upwinding function for the interface flux computation. If NULL,
      DefaultUpwinding() will be used. */
  int(*UpwindFunction)(double*,double*,double*,double*,
                       double*,double*,int,void*,double)
)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;

  double *gpu_uC     = NULL;
  double *gpu_uL     = solver->uL;
  double *gpu_uR     = solver->uR;
  double *gpu_fluxL  = solver->fL;
  double *gpu_fluxR  = solver->fR;

  /*
    precalculate the non-linear interpolation coefficients if required
    else reuse the weights previously calculated
  */
  if (LimFlag) solver->SetInterpLimiterVar(gpu_fluxC,gpu_u,gpu_x,dir,solver,mpi);

  /* if defined, calculate the modified u-function to be used for upwinding
     e.g.: used in well-balanced schemes for Euler/Navier-Stokes with gravity
     otherwise, just copy u to uC */
  if (solver->UFunction) {
    gpu_uC = solver->uC;
    solver->UFunction(gpu_uC,gpu_u,dir,solver,mpi,t);
  } else gpu_uC = gpu_u;

  /* Interpolation -> to calculate left and right-biased interface flux and state variable*/
  solver->InterpolateInterfacesHyp(gpu_uL   ,gpu_uC   ,gpu_u,gpu_x, 1,dir,solver,mpi,1);
  solver->InterpolateInterfacesHyp(gpu_uR   ,gpu_uC   ,gpu_u,gpu_x,-1,dir,solver,mpi,1);
  solver->InterpolateInterfacesHyp(gpu_fluxL,gpu_fluxC,gpu_u,gpu_x, 1,dir,solver,mpi,0);
  solver->InterpolateInterfacesHyp(gpu_fluxR,gpu_fluxC,gpu_u,gpu_x,-1,dir,solver,mpi,0);

  /* Upwind -> to calculate the final interface flux */
  if (UpwindFunction) { UpwindFunction   (gpu_fluxI,gpu_fluxL,gpu_fluxR,gpu_uL,gpu_uR,gpu_u,dir,solver,t); }
  else                { DefaultUpwinding (gpu_fluxI,gpu_fluxL,gpu_fluxR,NULL,NULL,NULL,dir,solver,t); }

  return(0);
}

#else

/*! 3D Kernel for gpuHyperbolicFunction() */
__global__
void HyperbolicFunction_dim3_kernel(
    int npoints_grid,
    int npoints_local_wghosts,
    int npoints_dim_interface,
    int d,
    int nvars,
    int ghosts,
    int offset,
    const int    * __restrict__ dim,
    const double * __restrict__ FluxI,
    const double * __restrict__ dxinv,
    double       * __restrict__ hyp,
    double       * __restrict__ FluxI_p1,
    double       * __restrict__ FluxI_p2
)
{
    int tx = threadIdx.x + (blockDim.x * blockIdx.x);
    if (tx < npoints_grid) {
        int v, p, p1, p2, b;
        int dim_interface[3];
        int index[3], index1[3], index2[3];

        _ArrayIndexnD_(3,tx,dim,index,0);
        _ArrayCopy1D_(dim,dim_interface,3);
        dim_interface[d] += 1;

        _ArrayCopy1D_(index,index1,3);
        _ArrayCopy1D_(index,index2,3); index2[d]++;
        _ArrayIndex1D_(3,dim          ,index ,ghosts,p);
        _ArrayIndex1D_(3,dim_interface,index1,0     ,p1);
        _ArrayIndex1D_(3,dim_interface,index2,0     ,p2);
        for (v=0; v<nvars; v++) {
          hyp[p+v*npoints_local_wghosts] += dxinv[offset+ghosts+index[d]]
                                      * (FluxI[p2+v*npoints_dim_interface]-FluxI[p1+v*npoints_dim_interface]);
        }
        if (index[d] == 0 || index[d] == dim[d]-1) {
          if (d == 0) {
            b = index[1] + index[2]*dim[1];
          } else if (d == 1) {
            b = index[0] + index[2]*dim[0];
          } else {
            b = index[0] + index[1]*dim[0];
          }
          if (index[d] == 0) {
            for (v=0; v<nvars; v++) FluxI_p1[b*nvars+v] = FluxI[p1+v*npoints_dim_interface];
          } else {
            for (v=0; v<nvars; v++) FluxI_p2[b*nvars+v] = FluxI[p2+v*npoints_dim_interface];
          }
        }
    }
    return;
}


/*! This function computes the hyperbolic term of the governing equations, expressed as follows:-
    \f{equation}{
      {\bf F} \left({\bf u}\right) = \sum_{d=0}^{D-1} \frac {\partial {\bf f}_d\left({\bf u}\right)} {\partial x_d}
    \f}
    using a conservative finite-difference discretization is space:
    \f{equation}{
      \hat{\bf F} \left({\bf u}\right) = \sum_{d=0}^{D-1} \frac {1} {\Delta x_d} \left[ \hat{\bf f}_{d,j+1/2} - \hat{\bf f}_{d,j-1/2} \right]
    \f}
    where \f$d\f$ denotes the spatial dimension, \f$D\f$ denotes the total number of spatial dimensions, the hat denotes
    the discretized quantity, and \f$j\f$ is the grid coordinate along dimension \f$d\f$.
    The approximation to the flux function \f${\bf f}_d\f$ at the interfaces \f$j\pm1/2\f$, denoted by \f$\hat{\bf f}_{d,j\pm 1/2}\f$,
    are computed using the function ReconstructHyperbolic().
*/
extern "C" int gpuHyperbolicFunction(
  double  *hyp, /*!< Array to hold the computed hyperbolic term (shares the same layout as u */
  double  *u,   /*!< Solution array */
  void    *s,      /*!< Solver object of type #HyPar */
  void    *m,      /*!< MPI object of type #MPIVariables */
  double  t,       /*!< Current simulation time */
  int     LimFlag, /*!< Flag to indicate if the nonlinear coefficients for solution-dependent
                          interpolation method should be recomputed (see ReconstructHyperbolic() for
                          an explanation on why this is needed) */
  /*! Function pointer to the flux function for the hyperbolic term */
  int(*FluxFunction)(double*,double*,int,void*,double),
  /*! Function pointer to the upwinding function for the hyperbolic term */
  int(*UpwindFunction)( double*,double*,double*,double*,
                        double*,double*,int,void*,double)
)
{
    HyPar         *solver = (HyPar*)        s;
    MPIVariables  *mpi    = (MPIVariables*) m;
    int           d;
    double        *FluxI  = solver->fluxI; /* interface flux     */
    double        *FluxC  = solver->fluxC; /* cell centered flux */

    int     ndims     = solver->ndims;
    int     nvars     = solver->nvars;
    int     ghosts    = solver->ghosts;
    int     *dim      = solver->dim_local;
    int     size      = solver->npoints_local_wghosts;
    double  *x        = solver->gpu_x;
    double  *dxinv    = solver->gpu_dxinv;
    double  *FluxI_p1, *FluxI_p2;

    LimFlag = (LimFlag && solver->flag_nonlinearinterp && solver->SetInterpLimiterVar);

    gpuMemset(hyp, 0, size*nvars*sizeof(double));
    gpuMemset(solver->StageBoundaryBuffer, 0, solver->StageBoundaryBuffer_size*sizeof(double));
    gpuMemset(solver->StageBoundaryIntegral, 0, 2*ndims*nvars*sizeof(double));

    if (!FluxFunction) return(0); /* zero hyperbolic term */
    /*solver->count_hyp++;*/

    int npoints_grid = 1; for (d = 0; d < ndims; d++) npoints_grid *= dim[d];
    int nblocks = (npoints_grid - 1) / GPU_THREADS_PER_BLOCK + 1;
    int offset = 0;

#if defined(GPU_STAT)
    cudaEvent_t start, stop;
    float milliseconds = 0;

    checkCuda( cudaEventCreate(&start) );
    checkCuda( cudaEventCreate(&stop) );
#endif

    for (d = 0; d < ndims; d++) {
        int npoints_dim_interface = 1; for (int i=0; i<ndims; i++) npoints_dim_interface *= (i==d) ? (dim[i]+1) : dim[i];
        /* evaluate cell-centered flux */
        FluxFunction(FluxC,u,d,solver,t);
        /* compute interface fluxes */
        ReconstructHyperbolic(FluxI,FluxC,u,x+offset,d,solver,mpi,t,LimFlag,UpwindFunction);

        FluxI_p1 = solver->StageBoundaryBuffer+(solver->gpu_npoints_boundary_offset[d]*nvars);
        FluxI_p2 = solver->StageBoundaryBuffer+(solver->gpu_npoints_boundary_offset[d]+solver->gpu_npoints_boundary[d])*nvars;

#if defined(GPU_STAT)
        checkCuda( cudaEventRecord(start, 0) );
#endif
        if (ndims == 3) {
          HyperbolicFunction_dim3_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
              npoints_grid, size, npoints_dim_interface, d, nvars, ghosts, offset,
              solver->gpu_dim_local, FluxI, dxinv,
              hyp, FluxI_p1, FluxI_p2
          );
        } else {
          fprintf(stderr,"gpuHyperbolicFunction(): Not implemented for ndims = %d!\n", ndims);
          exit(1);
        }
        cudaDeviceSynchronize();

#if defined(GPU_STAT)
        checkCuda( cudaEventRecord(stop, 0) );
        checkCuda( cudaEventSynchronize(stop) );
        checkCuda( cudaEventElapsedTime(&milliseconds, start, stop) );
        printf("%-50s GPU time = %.6f dir = %d bandwidth (GB/s) = %.2f\n",
                "HyperbolicFunction", milliseconds*1e-3, d, (1e-6*2*npoints_grid*nvars*sizeof(double))/milliseconds);

        checkCuda( cudaEventRecord(start, 0) );
#endif

        StageBoundaryIntegral_kernel<GPU_THREADS_PER_BLOCK><<<1, GPU_THREADS_PER_BLOCK, GPU_THREADS_PER_BLOCK*nvars*sizeof(double)>>>(
          solver->gpu_npoints_boundary[d], nvars, -1, FluxI_p1,
          solver->StageBoundaryIntegral + 2*d*nvars
        );
        StageBoundaryIntegral_kernel<GPU_THREADS_PER_BLOCK><<<1, GPU_THREADS_PER_BLOCK, GPU_THREADS_PER_BLOCK*nvars*sizeof(double)>>>(
          solver->gpu_npoints_boundary[d], nvars, 1, FluxI_p2,
          solver->StageBoundaryIntegral + (2*d+1)*nvars
        );
        cudaDeviceSynchronize();

#if defined(GPU_STAT)
        checkCuda( cudaEventRecord(stop, 0) );
        checkCuda( cudaEventSynchronize(stop) );
        checkCuda( cudaEventElapsedTime(&milliseconds, start, stop) );
        printf("%-50s GPU time = %.6f dir = %d bandwidth (GB/s) = %.2f\n",
                "StageBoundaryIntegral", milliseconds*1e-3, d, (1e-6*2*solver->gpu_npoints_boundary[d]*nvars*sizeof(double))/milliseconds);
#endif
        offset += dim[d] + 2*ghosts;
    }

    if (solver->flag_ib) gpuArrayBlockMultiply(hyp, solver->gpu_iblank, size, nvars);

#if defined(GPU_STAT)
    checkCuda(cudaEventDestroy(start));
    checkCuda(cudaEventDestroy(stop));
#endif

    return(0);
}

/*! This function computes the numerical flux \f$\hat{\bf f}_{j+1/2}\f$ at the interface from the cell-centered
    flux function \f${\bf f}_j\f$. This happens in two steps:-

    \b Interpolation: High-order accurate approximations to the flux at the interface \f$j+1/2\f$ are computed from
    the cell-centered flux with left- and right-biased interpolation methods. This is done by the
    #HyPar::InterpolateInterfacesHyp function. This can be expressed as follows:
    \f{align}{
      \hat{\bf f}^L_{j+1/2} &= \mathcal{I}\left({\bf f}_j,+1\right), \\
      \hat{\bf f}^R_{j+1/2} &= \mathcal{I}\left({\bf f}_j,-1\right),
    \f}
    where the \f$\pm 1\f$ indicates the interpolation bias, and \f$\mathcal{I}\f$ is the interpolation operator
    pointed to by #HyPar::InterpolateInterfacesHyp (see \b src/InterpolationFunctions for all the available operators).
    The interface values of the solution are similarly computed:
    \f{align}{
      \hat{\bf u}^L_{j+1/2} &= \mathcal{I}\left({\bf u}_j,+1\right), \\
      \hat{\bf u}^R_{j+1/2} &= \mathcal{I}\left({\bf u}_j,-1\right),
    \f}
    The specific choice of \f$\mathcal{I}\f$ is set based on #HyPar::spatial_scheme_hyp.

    \b Upwinding: The final flux at the interface is computed as
    \f{equation}{
      \hat{\bf f}_{j+1/2} = \mathcal{U}\left( \hat{\bf f}^L_{j+1/2}, \hat{\bf f}^R_{j+1/2}, \hat{\bf u}^L_{j+1/2}, \hat{\bf u}^R_{j+1/2} \right),
    \f}
    where \f$\mathcal{U}\f$ denotes the upwinding function UpwindFunction() passed as an argument (if NULL, DefaultUpwinding() is used). The
    upwinding function is specified by the physical model.

    \b Note:
    Solution-dependent, nonlinear interpolation methods (such as WENO, CRWENO) are implemented in a way that separates the calculation of the
    nonlinear interpolation weights (based on, say, the smoothness of the flux function), and the actual evaluation of the interpolant, into
    different functions. This allows the flexibility to choose if and when the nonlinear coefficients are evaluated (or previously computed
    values are reused). Some possible scenarios are:
    + For explicit time integration, they are computed every time the hyperbolic flux term is being computed.
    + For implicit time integration, consistency or linearization may require that they be computed and "frozen"
      at the beginning of a time step or stage.

    The argument \b LimFlag controls this behavior:
    + LimFlag = 1 means recompute the nonlinear coefficients.
    + LimFlag = 0 means reuse the the previously computed coefficients.
*/
int ReconstructHyperbolic(
  double  *fluxI,     /*!< Array to hold the computed interface fluxes. This array does not
                           have ghost points. The dimensions are the same as those of u without
                           ghost points in all dimensions, except along dir, where it is one more */
  double  *fluxC,     /*!< Array of the flux function computed at the cell centers
                           (same layout as u) */
  double  *u,         /*!< Solution array */
  double  *x,         /*!< Array of spatial coordinates */
  int     dir,        /*!< Spatial dimension along which to reconstruct the interface fluxes */
  void    *s,         /*!< Solver object of type #HyPar */
  void    *m,         /*!< MPI object of type #MPIVariables */
  double  t,          /*!< Current solution time */
  int     LimFlag,    /*!< Flag to indicate if the nonlinear coefficients for solution-dependent
                           interpolation method should be recomputed */
  /*! Function pointer to the upwinding function for the interface flux computation. If NULL,
      DefaultUpwinding() will be used. */
  int(*UpwindFunction)(double*,double*,double*,double*,
                       double*,double*,int,void*,double)
)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;

  double *uC     = NULL;
  double *uL     = solver->uL;
  double *uR     = solver->uR;
  double *fluxL  = solver->fL;
  double *fluxR  = solver->fR;

  /*
    precalculate the non-linear interpolation coefficients if required
    else reuse the weights previously calculated
  */

  if (LimFlag) solver->SetInterpLimiterVar(fluxC,u,x,dir,solver,mpi);

  /* if defined, calculate the modified u-function to be used for upwinding
     e.g.: used in well-balanced schemes for Euler/Navier-Stokes with gravity
     otherwise, just copy u to uC */
  if (solver->UFunction) {
    uC = solver->uC;
    solver->UFunction(uC,u,dir,solver,mpi,t);
  } else uC = u;


  /* Interpolation -> to calculate left and right-biased interface flux and state variable*/
  solver->InterpolateInterfacesHyp(uL   ,uC   ,u,x, 1,dir,solver,mpi,1);
  solver->InterpolateInterfacesHyp(uR   ,uC   ,u,x,-1,dir,solver,mpi,1);
  solver->InterpolateInterfacesHyp(fluxL,fluxC,u,x, 1,dir,solver,mpi,0);
  solver->InterpolateInterfacesHyp(fluxR,fluxC,u,x,-1,dir,solver,mpi,0);

  /* Upwind -> to calculate the final interface flux */
  if (UpwindFunction) { UpwindFunction   (fluxI,fluxL,fluxR,uL,uR,u,dir,solver,t); }
  else                { DefaultUpwinding (fluxI,fluxL,fluxR,NULL,NULL,NULL,dir,solver,t); }

  return(0);
}

#endif

/*! If no upwinding scheme is specified, this function defines the "upwind" flux as the
    arithmetic mean of the left- and right-biased fluxes. */
int DefaultUpwinding(
  double  *fI,  /*!< Computed upwind interface flux */
  double  *fL,  /*!< Left-biased reconstructed interface flux */
  double  *fR,  /*!< Right-biased reconstructed interface flux */
  double  *uL,  /*!< Left-biased reconstructed interface solution */
  double  *uR,  /*!< Right-biased reconstructed interface solution */
  double  *u,   /*!< Cell-centered solution */
  int     dir,  /*!< Spatial dimension */
  void    *s,   /*!< Solver object of type #HyPar */
  double  t     /*!< Current solution time */
)
{
  HyPar *solver = (HyPar*)    s;
  int   done;

  int *dim  = solver->dim_local;
  int ndims = solver->ndims;
  int nvars = solver->nvars;

  int bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  done = 0; int index_outer[ndims], index_inter[ndims];
  _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      int v; for (v=0; v<nvars; v++) fI[nvars*p+v] = 0.5 * (fL[nvars*p+v]+fR[nvars*p+v]);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
