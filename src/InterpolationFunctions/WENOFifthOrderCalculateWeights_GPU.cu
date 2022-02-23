/*! @file WENOFifthOrderCalculateWeights_GPU.cu
    @brief Functions to compute the nonlinear weights for WENO-type schemes
    @author Youngdae Kim
*/

#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <cuda_profiler_api.h>

#include <basic.h>
#include <basic_gpu.h>
#include <arrayfunctions.h>
#include <matmult_native.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

#include <arrayfunctions_gpu.h>


static int gpuWENOFifthOrderCalculateWeightsYC(double*,double*,double*,int,void*,void*);
static int gpuWENOFifthOrderCalculateWeightsM(double*,double*,double*,int,void*,void*);

/*! Compute the nonlinear weights for 5th order WENO-type schemes. This function is a wrapper that
    calls the appropriate function, depending on the type of WENO weights.
*/
extern "C" int gpuWENOFifthOrderCalculateWeights(
  double  *fC, /*!< Array of cell-centered values of the function \f${\bf f}\left({\bf u}\right)\f$ */
  double  *uC, /*!< Array of cell-centered values of the solution \f${\bf u}\f$ */
  double  *x,  /*!< Grid coordinates */
  int     dir, /*!< Spatial dimension along which to interpolation */
  void    *s,  /*!< Object of type #HyPar containing solver-related variables */
  void    *m   /*!< Object of type #MPIVariables containing MPI-related variables */
)
{
  HyPar           *solver = (HyPar*)          s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  MPIVariables    *mpi    = (MPIVariables*)   m;

  if (weno->yc)           return(gpuWENOFifthOrderCalculateWeightsYC (fC,uC,x,dir,solver,mpi));
  else if (weno->mapped)  return(gpuWENOFifthOrderCalculateWeightsM (fC,uC,x,dir,solver,mpi));
  else {
    printf("ERROR! WENO functions other than yc or mapped have not been implemented on GPUs.\n");
    return 1;
  }
}

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel for gpuWENOFifthOrderCalculateWeightsM() */
__global__
void WENOFifthOrderCalculateWeightsM_kernel(
    int ngrid_points,
    int ndims,
    int dir,
    int ghosts,
    int nvars,
    int weno_size,
    int offset,
    int stride,
    int is_crweno,
    int is_mpi_ip_zero,
    int is_mpi_ip_proc,
    double weno_eps,
    const int *dim,
    const double *fC,
    const double *uC,
    double *w1,
    double *w2,
    double *w3
)
{
    int p = threadIdx.x + (blockDim.x * blockIdx.x);
    if (p < ngrid_points) {
      const int    max_ndims = 3;
      const double thirteen_by_twelve = 13.0 / 12.0;
      const double one_fourth = 1.0 / 4.0;

      int    bounds_inter[max_ndims], indexC[max_ndims], indexI[max_ndims];
      int    qm1L, qm2L, qm3L, qp1L, qp2L, qm1R, qm2R, qm3R, qp1R, qp2R;
      double *ww1LF, *ww2LF, *ww3LF, *ww1RF, *ww2RF, *ww3RF;
      double *ww1LU, *ww2LU, *ww3LU, *ww1RU, *ww2RU, *ww3RU;

      ww1LF = w1 + offset;
      ww2LF = w2 + offset;
      ww3LF = w3 + offset;
      ww1RF = w1 + 2*weno_size + offset;
      ww2RF = w2 + 2*weno_size + offset;
      ww3RF = w3 + 2*weno_size + offset;
      ww1LU = w1 + weno_size + offset;
      ww2LU = w2 + weno_size + offset;
      ww3LU = w3 + weno_size + offset;
      ww1RU = w1 + 2*weno_size + weno_size + offset;
      ww2RU = w2 + 2*weno_size + weno_size + offset;
      ww3RU = w3 + 2*weno_size + weno_size + offset;

      _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
      _ArrayIndexnD_(ndims,p,bounds_inter,indexC,0);
      _ArrayCopy1D_(indexC,indexI,ndims);

      indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1L);
      qm3L = qm1L - 2*stride;
      qm2L = qm1L -   stride;
      qp1L = qm1L +   stride;
      qp2L = qm1L + 2*stride;

      indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1R);
      qm3R = qm1R + 2*stride;
      qm2R = qm1R +   stride;
      qp1R = qm1R -   stride;
      qp2R = qm1R - 2*stride;

      /* Defining stencil points */
      const double *m3LF, *m2LF, *m1LF, *p1LF, *p2LF;
      m3LF = (fC+qm3L*nvars);
      m2LF = (fC+qm2L*nvars);
      m1LF = (fC+qm1L*nvars);
      p1LF = (fC+qp1L*nvars);
      p2LF = (fC+qp2L*nvars);
      const double *m3RF, *m2RF, *m1RF, *p1RF, *p2RF;
      m3RF = (fC+qm3R*nvars);
      m2RF = (fC+qm2R*nvars);
      m1RF = (fC+qm1R*nvars);
      p1RF = (fC+qp1R*nvars);
      p2RF = (fC+qp2R*nvars);
      const double *m3LU, *m2LU, *m1LU, *p1LU, *p2LU;
      m3LU = (uC+qm3L*nvars);
      m2LU = (uC+qm2L*nvars);
      m1LU = (uC+qm1L*nvars);
      p1LU = (uC+qp1L*nvars);
      p2LU = (uC+qp2L*nvars);
      const double *m3RU, *m2RU, *m1RU, *p1RU, *p2RU;
      m3RU = (uC+qm3R*nvars);
      m2RU = (uC+qm2R*nvars);
      m1RU = (uC+qm1R*nvars);
      p1RU = (uC+qp1R*nvars);
      p2RU = (uC+qp2R*nvars);

      double c1, c2, c3;
      if (is_crweno) {
        if (   (is_mpi_ip_zero && (indexI[dir] == 0       ))
            || (is_mpi_ip_proc && (indexI[dir] == dim[dir])) ) {
          /* Use WENO5 at the physical boundaries */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        } else {
          /* CRWENO5 at the interior points */
          c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
          c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
          c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
        }
      } else {
        /* WENO5 and HCWENO5 */
        c1 = _WENO_OPTIMAL_WEIGHT_1_;
        c2 = _WENO_OPTIMAL_WEIGHT_2_;
        c3 = _WENO_OPTIMAL_WEIGHT_3_;
      }

      /* calculate WENO weights */
      _WENOWeights_v_M_((ww1LF+p*nvars),(ww2LF+p*nvars),(ww3LF+p*nvars),c1,c2,c3,m3LF,m2LF,m1LF,p1LF,p2LF,weno_eps,nvars);
      _WENOWeights_v_M_((ww1RF+p*nvars),(ww2RF+p*nvars),(ww3RF+p*nvars),c1,c2,c3,m3RF,m2RF,m1RF,p1RF,p2RF,weno_eps,nvars);
      _WENOWeights_v_M_((ww1LU+p*nvars),(ww2LU+p*nvars),(ww3LU+p*nvars),c1,c2,c3,m3LU,m2LU,m1LU,p1LU,p2LU,weno_eps,nvars);
      _WENOWeights_v_M_((ww1RU+p*nvars),(ww2RU+p*nvars),(ww3RU+p*nvars),c1,c2,c3,m3RU,m2RU,m1RU,p1RU,p2RU,weno_eps,nvars);
    }

    return;
}

/*!
  Computes the nonlinear weights for the 5th order component-wise WENO-type schemes using the mapped formulation of Henrick, Aslam & Powers:
  \f{eqnarray}{
    \omega_k &=& \frac {a_k} {\sum_{j=1}^3 a_j },\ a_k = \frac {\tilde{\omega}_k \left( c_k + c_k^2 - 3c_k\tilde{\omega}_k + \tilde{\omega}_k^2\right)} {c_k^2 + \tilde{\omega}_k\left(1-2c_k\right)}, \\
    \tilde{\omega}_k &=& \frac {\tilde{a}_k} {\sum_{j=1}^3 \tilde{a}_j },\ \tilde{a}_k = \frac {c_k} {\left(\beta_k+\epsilon\right)^p},\ k = 1,2,3,
  \f}
  where \f$c_k\f$ are the optimal weights, \f$p\f$ is hardcoded to \f$2\f$, and \f$\epsilon\f$ is an input parameter
  (#WENOParameters::eps) (typically \f$10^{-6}\f$). The smoothness indicators \f$\beta_k\f$ are given by:
  \f{eqnarray}{
    \beta_1 &=& \frac{13}{12} \left(f_{j-2}-2f_{j-1}+f_j\right)^2 + \frac{1}{4}\left(f_{j-2}-4f_{j-1}+3f_j\right)^2 \\
    \beta_2 &=& \frac{13}{12} \left(f_{j-1}-2f_j+f_{j+1}\right)^2 + \frac{1}{4}\left(f_{j-1}-f_{j+1}\right)^2 \\
    \beta_3 &=& \frac{13}{12} \left(f_j-2f_{j+1}+f_{j+2}\right)^2 + \frac{1}{4}\left(3f_j-4f_{j+1}+f_{j+2}\right)^2
  \f}

  \b Notes:
  + This function computes the weights for the entire grid (for interpolation along a given spatial dimension \a dir ).
    Thus, it loops over all grid lines along \a dir.
  + The weights are computed for all components, for both the left- and right-biased interpolations, and for both the
    flux function \f${\bf f}\left({\bf u}\right)\f$ and the solution \f$\bf u\f$. Thus, this approach of computing and
    storing the weights is quite memory-intensive.
  + The variable \a offset denotes the offset in the complete array from where the weights for interpolation along the
    current interpolation dimension (\a dir ) are stored.

 \b Reference:
 + Henrick, Aslam, Powers, Mapped weighted essentially non-oscillatory schemes: Achieving optimal order near critical
   points, J. Comput. Phys., 2005. http://dx.doi.org/10.1016/j.jcp.2005.01.023
*/
int gpuWENOFifthOrderCalculateWeightsM(
    double * __restrict__ fC, /*!< Array of cell-centered values of the function \f${\bf f}\left({\bf u}\right)\f$ */
    double * __restrict__ uC, /*!< Array of cell-centered values of the solution \f${\bf u}\f$ */
    double * __restrict__ x,  /*!< Grid coordinates */
    int dir,    /*!< Spatial dimension along which to interpolation */
    void *s,    /*!< Object of type #HyPar containing solver-related variables */
    void *m     /*!< Object of type #MPIVariables containing MPI-related variables */
)
{
    HyPar           *solver = (HyPar*)          s;
    WENOParameters  *weno   = (WENOParameters*) solver->interp;
    MPIVariables    *mpi    = (MPIVariables*)   m;

    int ghosts = solver->ghosts;
    int ndims  = solver->ndims;
    int nvars  = solver->nvars;
    int *dim   = solver->dim_local;
    int *stride= solver->stride_with_ghosts;

    /* calculate dimension offset */
    int offset = weno->offset[dir];
    int bounds_inter[ndims];
    _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] +=  1;
    int npoints_grid; _ArrayProduct1D_(bounds_inter,ndims,npoints_grid);
    int nblocks = (npoints_grid-1) / GPU_THREADS_PER_BLOCK + 1;

    int is_crweno      = (strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_) == 0) ? 1 : 0;
    int is_mpi_ip_zero = (mpi->ip[dir] == 0) ? 1 : 0;
    int is_mpi_ip_proc = (mpi->ip[dir] == mpi->iproc[dir]-1) ? 1 : 0;

#if defined(GPU_STAT)
    cudaEvent_t startEvent, stopEvent;
    float milliseconds = 0;

    checkCuda( cudaEventCreate(&startEvent) );
    checkCuda( cudaEventCreate(&stopEvent) );

    int weno_memory_accessed = 12*npoints_grid*nvars*sizeof(double);
    int fCuC_memory_accessed = 1;
    for (int d=0; d<ndims; d++) {
      if (d == dir) fCuC_memory_accessed *= (dim[d]+2*ghosts);
      else          fCuC_memory_accessed *= dim[d];
    }
    fCuC_memory_accessed *= nvars*sizeof(double);

    checkCuda( cudaEventRecord(startEvent, 0) );
#endif

    WENOFifthOrderCalculateWeightsM_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        npoints_grid, ndims, dir, ghosts, nvars, weno->size, offset, stride[dir],
        is_crweno, is_mpi_ip_zero, is_mpi_ip_proc, weno->eps,
        solver->gpu_dim_local, fC, uC, weno->w1, weno->w2, weno->w3
    );
    cudaDeviceSynchronize();

#if defined(GPU_STAT)
    checkCuda( cudaEventRecord(stopEvent, 0) );
    checkCuda( cudaEventSynchronize(stopEvent) );
    checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

    printf("%-50s GPU time (secs) = %.6f dir = %d bandwidth (GB/s) = %6.2f\n",
           "WENOFifthOrderCalculateWeightsM", milliseconds*1e-3, dir,
           (1e-6*(weno_memory_accessed+fCuC_memory_accessed)/milliseconds));

    checkCuda( cudaEventDestroy(startEvent) );
    checkCuda( cudaEventDestroy(stopEvent) );
#endif

    return 0;
}

/*! Kernel for gpuWENOFifthOrderCalculateWeightsYC() */
__global__
void WENOFifthOrderCalculateWeightsYC_kernel(
    int ngrid_points,
    int ndims,
    int dir,
    int ghosts,
    int nvars,
    int weno_size,
    int offset,
    int stride,
    int is_crweno,
    int is_mpi_ip_zero,
    int is_mpi_ip_proc,
    double weno_eps,
    const int *dim,
    const double *fC,
    const double *uC,
    double *w1,
    double *w2,
    double *w3
)
{
    int p = threadIdx.x + (blockDim.x * blockIdx.x);
    if (p < ngrid_points) {
      const int    max_ndims = 3;
      const double thirteen_by_twelve = 13.0 / 12.0;
      const double one_fourth = 1.0 / 4.0;

      int    bounds_inter[max_ndims], indexC[max_ndims], indexI[max_ndims];
      int    qm1L, qm2L, qm3L, qp1L, qp2L, qm1R, qm2R, qm3R, qp1R, qp2R;
      double *ww1LF, *ww2LF, *ww3LF, *ww1RF, *ww2RF, *ww3RF;
      double *ww1LU, *ww2LU, *ww3LU, *ww1RU, *ww2RU, *ww3RU;

      ww1LF = w1 + offset;
      ww2LF = w2 + offset;
      ww3LF = w3 + offset;
      ww1RF = w1 + 2*weno_size + offset;
      ww2RF = w2 + 2*weno_size + offset;
      ww3RF = w3 + 2*weno_size + offset;
      ww1LU = w1 + weno_size + offset;
      ww2LU = w2 + weno_size + offset;
      ww3LU = w3 + weno_size + offset;
      ww1RU = w1 + 2*weno_size + weno_size + offset;
      ww2RU = w2 + 2*weno_size + weno_size + offset;
      ww3RU = w3 + 2*weno_size + weno_size + offset;

      _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
      _ArrayIndexnD_(ndims,p,bounds_inter,indexC,0);
      _ArrayCopy1D_(indexC,indexI,ndims);

      indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1L);
      qm3L = qm1L - 2*stride;
      qm2L = qm1L -   stride;
      qp1L = qm1L +   stride;
      qp2L = qm1L + 2*stride;

      indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1R);
      qm3R = qm1R + 2*stride;
      qm2R = qm1R +   stride;
      qp1R = qm1R -   stride;
      qp2R = qm1R - 2*stride;

      /* Defining stencil points */
      const double *m3LF, *m2LF, *m1LF, *p1LF, *p2LF;
      m3LF = (fC+qm3L*nvars);
      m2LF = (fC+qm2L*nvars);
      m1LF = (fC+qm1L*nvars);
      p1LF = (fC+qp1L*nvars);
      p2LF = (fC+qp2L*nvars);
      const double *m3RF, *m2RF, *m1RF, *p1RF, *p2RF;
      m3RF = (fC+qm3R*nvars);
      m2RF = (fC+qm2R*nvars);
      m1RF = (fC+qm1R*nvars);
      p1RF = (fC+qp1R*nvars);
      p2RF = (fC+qp2R*nvars);
      const double *m3LU, *m2LU, *m1LU, *p1LU, *p2LU;
      m3LU = (uC+qm3L*nvars);
      m2LU = (uC+qm2L*nvars);
      m1LU = (uC+qm1L*nvars);
      p1LU = (uC+qp1L*nvars);
      p2LU = (uC+qp2L*nvars);
      const double *m3RU, *m2RU, *m1RU, *p1RU, *p2RU;
      m3RU = (uC+qm3R*nvars);
      m2RU = (uC+qm2R*nvars);
      m1RU = (uC+qm1R*nvars);
      p1RU = (uC+qp1R*nvars);
      p2RU = (uC+qp2R*nvars);

      double c1, c2, c3;
      if (is_crweno) {
        if (   (is_mpi_ip_zero && (indexI[dir] == 0       ))
            || (is_mpi_ip_proc && (indexI[dir] == dim[dir])) ) {
          /* Use WENO5 at the physical boundaries */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        } else {
          /* CRWENO5 at the interior points */
          c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
          c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
          c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
        }
      } else {
        /* WENO5 and HCWENO5 */
        c1 = _WENO_OPTIMAL_WEIGHT_1_;
        c2 = _WENO_OPTIMAL_WEIGHT_2_;
        c3 = _WENO_OPTIMAL_WEIGHT_3_;
      }

      /* calculate WENO weights */
      _WENOWeights_v_YC_((ww1LF+p*nvars),(ww2LF+p*nvars),(ww3LF+p*nvars),c1,c2,c3,m3LF,m2LF,m1LF,p1LF,p2LF,weno_eps,nvars);
      _WENOWeights_v_YC_((ww1RF+p*nvars),(ww2RF+p*nvars),(ww3RF+p*nvars),c1,c2,c3,m3RF,m2RF,m1RF,p1RF,p2RF,weno_eps,nvars);
      _WENOWeights_v_YC_((ww1LU+p*nvars),(ww2LU+p*nvars),(ww3LU+p*nvars),c1,c2,c3,m3LU,m2LU,m1LU,p1LU,p2LU,weno_eps,nvars);
      _WENOWeights_v_YC_((ww1RU+p*nvars),(ww2RU+p*nvars),(ww3RU+p*nvars),c1,c2,c3,m3RU,m2RU,m1RU,p1RU,p2RU,weno_eps,nvars);
    }

    return;
}

/*!
  Computes the nonlinear weights for the 5th order component-wise WENO-type schemes using the ESWENO formulation of
  Yamaleev & Carpenter. Note that only the formulation for the nonlinear weights is adopted and implemented here, not
  the ESWENO scheme as a whole.
  \f{equation}{
    \omega_k = \frac {a_k} {\sum_{j=1}^3 a_j },\ a_k = c_k \left( 1 + \frac{\tau_5}{\beta_k+\epsilon} \right)^p,\ k = 1,2,3,
  \f}
  where \f$c_k\f$ are the optimal weights, \f$p\f$ is hardcoded to \f$2\f$, and \f$\epsilon\f$ is an input parameter
  (#WENOParameters::eps) (typically \f$10^{-6}\f$). The smoothness indicators \f$\beta_k\f$ are given by:
  \f{eqnarray}{
    \beta_1 &=& \frac{13}{12} \left(f_{j-2}-2f_{j-1}+f_j\right)^2 + \frac{1}{4}\left(f_{j-2}-4f_{j-1}+3f_j\right)^2 \\
    \beta_2 &=& \frac{13}{12} \left(f_{j-1}-2f_j+f_{j+1}\right)^2 + \frac{1}{4}\left(f_{j-1}-f_{j+1}\right)^2 \\
    \beta_3 &=& \frac{13}{12} \left(f_j-2f_{j+1}+f_{j+2}\right)^2 + \frac{1}{4}\left(3f_j-4f_{j+1}+f_{j+2}\right)^2,
  \f}
  and \f$\tau_5 = \left( f_{j-2}-4f_{j-1}+6f_j-4f_{j+1}+f_{j+2} \right)^2\f$.

  \b Notes:
  + This function computes the weights for the entire grid (for interpolation along a given spatial dimension \a dir ).
    Thus, it loops over all grid lines along \a dir.
  + The weights are computed for all components, for both the left- and right-biased interpolations, and for both the
    flux function \f${\bf f}\left({\bf u}\right)\f$ and the solution \f$\bf u\f$. Thus, this approach of computing and
    storing the weights is quite memory-intensive.
  + The variable \a offset denotes the offset in the complete array from where the weights for interpolation along the
    current interpolation dimension (\a dir ) are stored.

  \b Reference:
     + Yamaleev, Carpenter, A systematic methodology for constructing high-order energy stable WENO schemes,
       J. Comput. Phys., 2009. http://dx.doi.org/10.1016/j.jcp.2009.03.002
*/
int gpuWENOFifthOrderCalculateWeightsYC(
    double * __restrict__ fC, /*!< Array of cell-centered values of the function \f${\bf f}\left({\bf u}\right)\f$ */
    double * __restrict__ uC, /*!< Array of cell-centered values of the solution \f${\bf u}\f$ */
    double * __restrict__ x,  /*!< Grid coordinates */
    int dir,    /*!< Spatial dimension along which to interpolation */
    void *s,    /*!< Object of type #HyPar containing solver-related variables */
    void *m     /*!< Object of type #MPIVariables containing MPI-related variables */
)
{
    HyPar           *solver = (HyPar*)          s;
    WENOParameters  *weno   = (WENOParameters*) solver->interp;
    MPIVariables    *mpi    = (MPIVariables*)   m;

    int ghosts = solver->ghosts;
    int ndims  = solver->ndims;
    int nvars  = solver->nvars;
    int *dim   = solver->dim_local;
    int *stride= solver->stride_with_ghosts;

    /* calculate dimension offset */
    int offset = weno->offset[dir];
    int bounds_inter[ndims];
    _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] +=  1;
    int ngrid_points; _ArrayProduct1D_(bounds_inter,ndims,ngrid_points);
    //int nblocks = (ngrid_points - 1) / GPU_THREADS_PER_BLOCK + 1;
    int nblocks = (ngrid_points - 1) / 256 + 1;

    int is_crweno      = (strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_) == 0) ? 1 : 0;
    int is_mpi_ip_zero = (mpi->ip[dir] == 0) ? 1 : 0;
    int is_mpi_ip_proc = (mpi->ip[dir] == mpi->iproc[dir]-1) ? 1 : 0;

    cudaEvent_t startEvent, stopEvent;
    float milliseconds = 0;

    checkCuda( cudaEventCreate(&startEvent) );
    checkCuda( cudaEventCreate(&stopEvent) );

    //cudaProfilerStart();
    checkCuda( cudaEventRecord(startEvent, 0) );
    WENOFifthOrderCalculateWeightsYC_kernel<<<nblocks, 256>>>(
        ngrid_points, ndims, dir, ghosts, nvars, weno->size, offset, stride[dir],
        is_crweno, is_mpi_ip_zero, is_mpi_ip_proc, weno->eps,
        solver->gpu_dim_local, fC, uC, weno->w1, weno->w2, weno->w3
    );
    checkCuda( cudaEventRecord(stopEvent, 0) );
    checkCuda( cudaEventSynchronize(stopEvent) );
    checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );
    //cudaProfilerStop();

    printf("WENOFifthOrderCalculateWeightsYC_GPU:\n");
    printf("  GPU time = %.6f dir = %d bandwidth (GB/s) = %.2f weno->size = %d weno->offset[1] = %d\n",
            milliseconds*1e-3, dir,
            (1e-6*((12*weno->offset[1]+2*(nvars*(dim[0]+2*ghosts)*dim[1]))*sizeof(double)))/milliseconds,
            weno->size, weno->offset[1]);

    checkCuda( cudaEventDestroy(startEvent) );
    checkCuda( cudaEventDestroy(stopEvent) );

    return (0);
}

#else

/*! Kernel for gpuWENOFifthOrderCalculateWeightsM() */
__global__
void WENOFifthOrderCalculateWeightsM_kernel(
    int npoints_grid,
    int npoints_local_wghosts,
    int ndims,
    int dir,
    int ghosts,
    int nvars,
    int weno_size,
    int offset,
    int stride,
    int is_crweno,
    int is_mpi_ip_zero,
    int is_mpi_ip_proc,
    double weno_eps,
    const int *dim,
    const double *fC,
    const double *uC,
    double *w1,
    double *w2,
    double *w3
)
{
    int p = threadIdx.x + (blockDim.x * blockIdx.x);
    if (p < npoints_grid) {
      const int    max_ndims = 3;
      const double thirteen_by_twelve = 13.0 / 12.0;
      const double one_fourth = 1.0 / 4.0;

      int    bounds_inter[max_ndims], indexC[max_ndims], indexI[max_ndims];
      int    qm1L, qm2L, qm3L, qp1L, qp2L, qm1R, qm2R, qm3R, qp1R, qp2R;
      double *ww1LF, *ww2LF, *ww3LF, *ww1RF, *ww2RF, *ww3RF;
      double *ww1LU, *ww2LU, *ww3LU, *ww1RU, *ww2RU, *ww3RU;

      ww1LF = w1 + offset;
      ww2LF = w2 + offset;
      ww3LF = w3 + offset;
      ww1RF = w1 + 2*weno_size + offset;
      ww2RF = w2 + 2*weno_size + offset;
      ww3RF = w3 + 2*weno_size + offset;
      ww1LU = w1 + weno_size + offset;
      ww2LU = w2 + weno_size + offset;
      ww3LU = w3 + weno_size + offset;
      ww1RU = w1 + 2*weno_size + weno_size + offset;
      ww2RU = w2 + 2*weno_size + weno_size + offset;
      ww3RU = w3 + 2*weno_size + weno_size + offset;

      _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
      _ArrayIndexnD_(ndims,p,bounds_inter,indexC,0);
      _ArrayCopy1D_(indexC,indexI,ndims);

      indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1L);
      qm3L = qm1L - 2*stride;
      qm2L = qm1L -   stride;
      qp1L = qm1L +   stride;
      qp2L = qm1L + 2*stride;

      indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1R);
      qm3R = qm1R + 2*stride;
      qm2R = qm1R +   stride;
      qp1R = qm1R -   stride;
      qp2R = qm1R - 2*stride;

      for (int j=0; j<nvars; j++) {
        /* Defining stencil points */
        const double *m3LF, *m2LF, *m1LF, *p1LF, *p2LF;
        m3LF = (fC+qm3L);
        m2LF = (fC+qm2L);
        m1LF = (fC+qm1L);
        p1LF = (fC+qp1L);
        p2LF = (fC+qp2L);
        const double *m3RF, *m2RF, *m1RF, *p1RF, *p2RF;
        m3RF = (fC+qm3R);
        m2RF = (fC+qm2R);
        m1RF = (fC+qm1R);
        p1RF = (fC+qp1R);
        p2RF = (fC+qp2R);
        const double *m3LU, *m2LU, *m1LU, *p1LU, *p2LU;
        m3LU = (uC+qm3L);
        m2LU = (uC+qm2L);
        m1LU = (uC+qm1L);
        p1LU = (uC+qp1L);
        p2LU = (uC+qp2L);
        const double *m3RU, *m2RU, *m1RU, *p1RU, *p2RU;
        m3RU = (uC+qm3R);
        m2RU = (uC+qm2R);
        m1RU = (uC+qm1R);
        p1RU = (uC+qp1R);
        p2RU = (uC+qp2R);

        qm3L += npoints_local_wghosts;
        qm2L += npoints_local_wghosts;
        qm1L += npoints_local_wghosts;
        qp1L += npoints_local_wghosts;
        qp2L += npoints_local_wghosts;

        qm3R += npoints_local_wghosts;
        qm2R += npoints_local_wghosts;
        qm1R += npoints_local_wghosts;
        qp1R += npoints_local_wghosts;
        qp2R += npoints_local_wghosts;

        double c1, c2, c3;
        if (is_crweno) {
          if (   (is_mpi_ip_zero && (indexI[dir] == 0       ))
              || (is_mpi_ip_proc && (indexI[dir] == dim[dir])) ) {
            /* Use WENO5 at the physical boundaries */
            c1 = _WENO_OPTIMAL_WEIGHT_1_;
            c2 = _WENO_OPTIMAL_WEIGHT_2_;
            c3 = _WENO_OPTIMAL_WEIGHT_3_;
          } else {
            /* CRWENO5 at the interior points */
            c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
            c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
            c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
          }
        } else {
          /* WENO5 and HCWENO5 */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        }

        /* calculate WENO weights */
        _WENOWeights_v_M_Scalar_((ww1LF+p),(ww2LF+p),(ww3LF+p),c1,c2,c3,m3LF,m2LF,m1LF,p1LF,p2LF,weno_eps,0);
        _WENOWeights_v_M_Scalar_((ww1RF+p),(ww2RF+p),(ww3RF+p),c1,c2,c3,m3RF,m2RF,m1RF,p1RF,p2RF,weno_eps,0);
        _WENOWeights_v_M_Scalar_((ww1LU+p),(ww2LU+p),(ww3LU+p),c1,c2,c3,m3LU,m2LU,m1LU,p1LU,p2LU,weno_eps,0);
        _WENOWeights_v_M_Scalar_((ww1RU+p),(ww2RU+p),(ww3RU+p),c1,c2,c3,m3RU,m2RU,m1RU,p1RU,p2RU,weno_eps,0);
        p += npoints_grid;
      }
    }

    return;
}


/*!
  Computes the nonlinear weights for the 5th order component-wise WENO-type schemes using the mapped formulation of Henrick, Aslam & Powers:
  \f{eqnarray}{
    \omega_k &=& \frac {a_k} {\sum_{j=1}^3 a_j },\ a_k = \frac {\tilde{\omega}_k \left( c_k + c_k^2 - 3c_k\tilde{\omega}_k + \tilde{\omega}_k^2\right)} {c_k^2 + \tilde{\omega}_k\left(1-2c_k\right)}, \\
    \tilde{\omega}_k &=& \frac {\tilde{a}_k} {\sum_{j=1}^3 \tilde{a}_j },\ \tilde{a}_k = \frac {c_k} {\left(\beta_k+\epsilon\right)^p},\ k = 1,2,3,
  \f}
  where \f$c_k\f$ are the optimal weights, \f$p\f$ is hardcoded to \f$2\f$, and \f$\epsilon\f$ is an input parameter
  (#WENOParameters::eps) (typically \f$10^{-6}\f$). The smoothness indicators \f$\beta_k\f$ are given by:
  \f{eqnarray}{
    \beta_1 &=& \frac{13}{12} \left(f_{j-2}-2f_{j-1}+f_j\right)^2 + \frac{1}{4}\left(f_{j-2}-4f_{j-1}+3f_j\right)^2 \\
    \beta_2 &=& \frac{13}{12} \left(f_{j-1}-2f_j+f_{j+1}\right)^2 + \frac{1}{4}\left(f_{j-1}-f_{j+1}\right)^2 \\
    \beta_3 &=& \frac{13}{12} \left(f_j-2f_{j+1}+f_{j+2}\right)^2 + \frac{1}{4}\left(3f_j-4f_{j+1}+f_{j+2}\right)^2
  \f}

  \b Notes:
  + This function computes the weights for the entire grid (for interpolation along a given spatial dimension \a dir ).
    Thus, it loops over all grid lines along \a dir.
  + The weights are computed for all components, for both the left- and right-biased interpolations, and for both the
    flux function \f${\bf f}\left({\bf u}\right)\f$ and the solution \f$\bf u\f$. Thus, this approach of computing and
    storing the weights is quite memory-intensive.
  + The variable \a offset denotes the offset in the complete array from where the weights for interpolation along the
    current interpolation dimension (\a dir ) are stored.

 \b Reference:
 + Henrick, Aslam, Powers, Mapped weighted essentially non-oscillatory schemes: Achieving optimal order near critical
   points, J. Comput. Phys., 2005. http://dx.doi.org/10.1016/j.jcp.2005.01.023
*/
int gpuWENOFifthOrderCalculateWeightsM(
    double * __restrict__ fC, /*!< Array of cell-centered values of the function \f${\bf f}\left({\bf u}\right)\f$ */
    double * __restrict__ uC, /*!< Array of cell-centered values of the solution \f${\bf u}\f$ */
    double * __restrict__ x,  /*!< Grid coordinates */
    int dir,    /*!< Spatial dimension along which to interpolation */
    void *s,    /*!< Object of type #HyPar containing solver-related variables */
    void *m     /*!< Object of type #MPIVariables containing MPI-related variables */
)
{
    HyPar           *solver = (HyPar*)          s;
    WENOParameters  *weno   = (WENOParameters*) solver->interp;
    MPIVariables    *mpi    = (MPIVariables*)   m;

    int ghosts = solver->ghosts;
    int ndims  = solver->ndims;
    int nvars  = solver->nvars;
    int *dim   = solver->dim_local;
    int *stride= solver->stride_with_ghosts;

    /* calculate dimension offset */
    int offset = weno->offset[dir];
    int bounds_inter[ndims];
    _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] +=  1;
    int npoints_grid; _ArrayProduct1D_(bounds_inter,ndims,npoints_grid);
    int nblocks = (npoints_grid-1) / GPU_THREADS_PER_BLOCK + 1;

    int is_crweno      = (strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_) == 0) ? 1 : 0;
    int is_mpi_ip_zero = (mpi->ip[dir] == 0) ? 1 : 0;
    int is_mpi_ip_proc = (mpi->ip[dir] == mpi->iproc[dir]-1) ? 1 : 0;

#if defined(GPU_STAT)
    cudaEvent_t startEvent, stopEvent;
    float milliseconds = 0;

    checkCuda( cudaEventCreate(&startEvent) );
    checkCuda( cudaEventCreate(&stopEvent) );

    int weno_memory_accessed = 12*npoints_grid*nvars*sizeof(double);
    int fCuC_memory_accessed = 1;
    for (int d=0; d<ndims; d++) {
      if (d == dir) fCuC_memory_accessed *= (dim[d]+2*ghosts);
      else          fCuC_memory_accessed *= dim[d];
    }
    fCuC_memory_accessed *= nvars*sizeof(double);

    checkCuda( cudaEventRecord(startEvent, 0) );
#endif

    WENOFifthOrderCalculateWeightsM_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        npoints_grid, solver->npoints_local_wghosts, ndims, dir, ghosts, nvars, weno->size, offset, stride[dir],
        is_crweno, is_mpi_ip_zero, is_mpi_ip_proc, weno->eps,
        solver->gpu_dim_local, fC, uC, weno->w1, weno->w2, weno->w3
    );
    cudaDeviceSynchronize();

#if defined(GPU_STAT)
    checkCuda( cudaEventRecord(stopEvent, 0) );
    checkCuda( cudaEventSynchronize(stopEvent) );
    checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

    printf("%-50s GPU time (secs) = %.6f dir = %d bandwidth (GB/s) = %6.2f stride = %d\n",
           "WENOFifthOrderCalculateWeightsM2", milliseconds*1e-3, dir,
           (1e-6*(weno_memory_accessed+fCuC_memory_accessed)/milliseconds),
           stride[dir]);

    checkCuda( cudaEventDestroy(startEvent) );
    checkCuda( cudaEventDestroy(stopEvent) );
#endif

    return 0;
}

/*! Kernel for gpuWENOFifthOrderCalculateWeightsYC() */
__global__
void WENOFifthOrderCalculateWeightsYC_kernel(
    int npoints_grid,
    int npoints_local_wghosts,
    int ndims,
    int dir,
    int ghosts,
    int nvars,
    int weno_size,
    int offset,
    int stride,
    int is_crweno,
    int is_mpi_ip_zero,
    int is_mpi_ip_proc,
    double weno_eps,
    const int * __restrict__ dim,
    const double * __restrict__ fC,
    const double * __restrict__ uC,
    double * __restrict__ w1,
    double * __restrict__ w2,
    double * __restrict__ w3
)
{
    int p = threadIdx.x + (blockDim.x * blockIdx.x);
    if (p < npoints_grid) {
      const int    max_ndims = 3;
      const double thirteen_by_twelve = 13.0 / 12.0;
      const double one_fourth = 1.0 / 4.0;

      int    bounds_inter[max_ndims], indexC[max_ndims], indexI[max_ndims];
      int    qm1L, qm2L, qm3L, qp1L, qp2L, qm1R, qm2R, qm3R, qp1R, qp2R;
      double *ww1LF, *ww2LF, *ww3LF, *ww1RF, *ww2RF, *ww3RF;
      double *ww1LU, *ww2LU, *ww3LU, *ww1RU, *ww2RU, *ww3RU;

      ww1LF = w1 + offset;
      ww2LF = w2 + offset;
      ww3LF = w3 + offset;
      ww1RF = w1 + 2*weno_size + offset;
      ww2RF = w2 + 2*weno_size + offset;
      ww3RF = w3 + 2*weno_size + offset;
      ww1LU = w1 + weno_size + offset;
      ww2LU = w2 + weno_size + offset;
      ww3LU = w3 + weno_size + offset;
      ww1RU = w1 + 2*weno_size + weno_size + offset;
      ww2RU = w2 + 2*weno_size + weno_size + offset;
      ww3RU = w3 + 2*weno_size + weno_size + offset;

      _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
      _ArrayIndexnD_(ndims,p,bounds_inter,indexC,0);
      _ArrayCopy1D_(indexC,indexI,ndims);

      indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1L);
      qm3L = qm1L - 2*stride;
      qm2L = qm1L -   stride;
      qp1L = qm1L +   stride;
      qp2L = qm1L + 2*stride;

      indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1R);
      qm3R = qm1R + 2*stride;
      qm2R = qm1R +   stride;
      qp1R = qm1R -   stride;
      qp2R = qm1R - 2*stride;

      for (int j=0; j <nvars; j++) {
        /* Defining stencil points */
        const double *m3LF, *m2LF, *m1LF, *p1LF, *p2LF;
        m3LF = (fC+qm3L);
        m2LF = (fC+qm2L);
        m1LF = (fC+qm1L);
        p1LF = (fC+qp1L);
        p2LF = (fC+qp2L);
        const double *m3RF, *m2RF, *m1RF, *p1RF, *p2RF;
        m3RF = (fC+qm3R);
        m2RF = (fC+qm2R);
        m1RF = (fC+qm1R);
        p1RF = (fC+qp1R);
        p2RF = (fC+qp2R);
        const double *m3LU, *m2LU, *m1LU, *p1LU, *p2LU;
        m3LU = (uC+qm3L);
        m2LU = (uC+qm2L);
        m1LU = (uC+qm1L);
        p1LU = (uC+qp1L);
        p2LU = (uC+qp2L);
        const double *m3RU, *m2RU, *m1RU, *p1RU, *p2RU;
        m3RU = (uC+qm3R);
        m2RU = (uC+qm2R);
        m1RU = (uC+qm1R);
        p1RU = (uC+qp1R);
        p2RU = (uC+qp2R);

        qm3L += npoints_local_wghosts;
        qm2L += npoints_local_wghosts;
        qm1L += npoints_local_wghosts;
        qp1L += npoints_local_wghosts;
        qp2L += npoints_local_wghosts;

        qm3R += npoints_local_wghosts;
        qm2R += npoints_local_wghosts;
        qm1R += npoints_local_wghosts;
        qp1R += npoints_local_wghosts;
        qp2R += npoints_local_wghosts;

        double c1, c2, c3;
        if (is_crweno) {
          if (   (is_mpi_ip_zero && (indexI[dir] == 0       ))
              || (is_mpi_ip_proc && (indexI[dir] == dim[dir])) ) {
            /* Use WENO5 at the physical boundaries */
            c1 = _WENO_OPTIMAL_WEIGHT_1_;
            c2 = _WENO_OPTIMAL_WEIGHT_2_;
            c3 = _WENO_OPTIMAL_WEIGHT_3_;
          } else {
            /* CRWENO5 at the interior points */
            c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
            c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
            c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
          }
        } else {
          /* WENO5 and HCWENO5 */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        }

        /* calculate WENO weights */
        _WENOWeights_v_YC_Scalar_((ww1LF+p),(ww2LF+p),(ww3LF+p),c1,c2,c3,m3LF,m2LF,m1LF,p1LF,p2LF,weno_eps,0);
        _WENOWeights_v_YC_Scalar_((ww1RF+p),(ww2RF+p),(ww3RF+p),c1,c2,c3,m3RF,m2RF,m1RF,p1RF,p2RF,weno_eps,0);
        _WENOWeights_v_YC_Scalar_((ww1LU+p),(ww2LU+p),(ww3LU+p),c1,c2,c3,m3LU,m2LU,m1LU,p1LU,p2LU,weno_eps,0);
        _WENOWeights_v_YC_Scalar_((ww1RU+p),(ww2RU+p),(ww3RU+p),c1,c2,c3,m3RU,m2RU,m1RU,p1RU,p2RU,weno_eps,0);
        p += npoints_grid;
      }
    }

    return;
}

/*!
  Computes the nonlinear weights for the 5th order component-wise WENO-type schemes using the ESWENO formulation of
  Yamaleev & Carpenter. Note that only the formulation for the nonlinear weights is adopted and implemented here, not
  the ESWENO scheme as a whole.
  \f{equation}{
    \omega_k = \frac {a_k} {\sum_{j=1}^3 a_j },\ a_k = c_k \left( 1 + \frac{\tau_5}{\beta_k+\epsilon} \right)^p,\ k = 1,2,3,
  \f}
  where \f$c_k\f$ are the optimal weights, \f$p\f$ is hardcoded to \f$2\f$, and \f$\epsilon\f$ is an input parameter
  (#WENOParameters::eps) (typically \f$10^{-6}\f$). The smoothness indicators \f$\beta_k\f$ are given by:
  \f{eqnarray}{
    \beta_1 &=& \frac{13}{12} \left(f_{j-2}-2f_{j-1}+f_j\right)^2 + \frac{1}{4}\left(f_{j-2}-4f_{j-1}+3f_j\right)^2 \\
    \beta_2 &=& \frac{13}{12} \left(f_{j-1}-2f_j+f_{j+1}\right)^2 + \frac{1}{4}\left(f_{j-1}-f_{j+1}\right)^2 \\
    \beta_3 &=& \frac{13}{12} \left(f_j-2f_{j+1}+f_{j+2}\right)^2 + \frac{1}{4}\left(3f_j-4f_{j+1}+f_{j+2}\right)^2,
  \f}
  and \f$\tau_5 = \left( f_{j-2}-4f_{j-1}+6f_j-4f_{j+1}+f_{j+2} \right)^2\f$.

  \b Notes:
  + This function computes the weights for the entire grid (for interpolation along a given spatial dimension \a dir ).
    Thus, it loops over all grid lines along \a dir.
  + The weights are computed for all components, for both the left- and right-biased interpolations, and for both the
    flux function \f${\bf f}\left({\bf u}\right)\f$ and the solution \f$\bf u\f$. Thus, this approach of computing and
    storing the weights is quite memory-intensive.
  + The variable \a offset denotes the offset in the complete array from where the weights for interpolation along the
    current interpolation dimension (\a dir ) are stored.

  \b Reference:
     + Yamaleev, Carpenter, A systematic methodology for constructing high-order energy stable WENO schemes,
       J. Comput. Phys., 2009. http://dx.doi.org/10.1016/j.jcp.2009.03.002
*/
int gpuWENOFifthOrderCalculateWeightsYC(
    double * __restrict__ fC, /*!< Array of cell-centered values of the function \f${\bf f}\left({\bf u}\right)\f$ */
    double * __restrict__ uC, /*!< Array of cell-centered values of the solution \f${\bf u}\f$ */
    double * __restrict__ x,  /*!< Grid coordinates */
    int dir,    /*!< Spatial dimension along which to interpolation */
    void *s,    /*!< Object of type #HyPar containing solver-related variables */
    void *m     /*!< Object of type #MPIVariables containing MPI-related variables */
)
{
    HyPar           *solver = (HyPar*)          s;
    WENOParameters  *weno   = (WENOParameters*) solver->interp;
    MPIVariables    *mpi    = (MPIVariables*)   m;

    int ghosts = solver->ghosts;
    int ndims  = solver->ndims;
    int nvars  = solver->nvars;
    int *dim   = solver->dim_local;
    int *stride= solver->stride_with_ghosts;

    /* calculate dimension offset */
    int offset = weno->offset[dir];
    int bounds_inter[ndims];
    _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] +=  1;
    int npoints_grid; _ArrayProduct1D_(bounds_inter,ndims,npoints_grid);
    int nblocks = (npoints_grid-1) / GPU_THREADS_PER_BLOCK + 1;

    int is_crweno      = (strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_) == 0) ? 1 : 0;
    int is_mpi_ip_zero = (mpi->ip[dir] == 0) ? 1 : 0;
    int is_mpi_ip_proc = (mpi->ip[dir] == mpi->iproc[dir]-1) ? 1 : 0;

#if defined(GPU_STAT)
    cudaEvent_t startEvent, stopEvent;
    float milliseconds = 0;

    checkCuda( cudaEventCreate(&startEvent) );
    checkCuda( cudaEventCreate(&stopEvent) );

    int weno_memory_accessed = 12*npoints_grid*nvars*sizeof(double);
    int fCuC_memory_accessed = 1;
    for (int d=0; d<ndims; d++) {
      if (d == dir) fCuC_memory_accessed *= (dim[d]+2*ghosts);
      else          fCuC_memory_accessed *= dim[d];
    }
    fCuC_memory_accessed *= nvars*sizeof(double);

    checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  WENOFifthOrderCalculateWeightsYC_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      npoints_grid, solver->npoints_local_wghosts, ndims, dir, ghosts, nvars, weno->size, offset, stride[dir],
      is_crweno, is_mpi_ip_zero, is_mpi_ip_proc, weno->eps,
      solver->gpu_dim_local, fC, uC, weno->w1, weno->w2, weno->w3
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
    checkCuda( cudaEventRecord(stopEvent, 0) );
    checkCuda( cudaEventSynchronize(stopEvent) );
    checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

    printf("%-50s GPU time (secs) = %.6f dir = %d bandwidth (GB/s) = %6.2f stride = %d\n",
            "WENOFifthOrderCalculateWeightsYC2", milliseconds*1e-3, dir,
            (1e-6*(weno_memory_accessed+fCuC_memory_accessed)/milliseconds),
            stride[dir]);

    checkCuda( cudaEventDestroy(startEvent) );
    checkCuda( cudaEventDestroy(stopEvent) );
#endif

    return (0);
}

#endif

