/*! @file Interp1PrimFifthOrderWENO_GPU.cu
 *  @brief WENO5 Scheme (Component-wise application to vectors).
 *  @author Youngdae Kim
*/

#include <basic.h>
#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

#undef  _MINIMUM_GHOSTS_
/*! \def _MINIMUM_GHOSTS_
 * Minimum number of ghost points required for this interpolation
 * method.
*/
#define _MINIMUM_GHOSTS_ 3

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel for gpuInterp1PrimFifthOrderWENO() */
__global__
void Interp1PrimFifthOrderWENO_kernel(
    int ngrid_points,
    int ndims,
    int dir,
    int ghosts,
    int nvars,
    int weno_size,
    int offset,
    int stride,
    int upw,
    int uflag,
    const int *dim,
    const double *fC,
    const double *w1,
    const double *w2,
    const double *w3,
    double *fI
)
{
    int p = threadIdx.x + (blockDim.x * blockIdx.x);
    if (p < ngrid_points) {
      int    bounds_inter[GPU_MAX_NDIMS], indexC[GPU_MAX_NDIMS], indexI[GPU_MAX_NDIMS];
      int qm1,qm2,qm3,qp1,qp2;
      const double *ww1, *ww2, *ww3;

      ww1 = w1 + (upw < 0 ? 2*weno_size : 0) + (uflag ? weno_size : 0) + offset;
      ww2 = w2 + (upw < 0 ? 2*weno_size : 0) + (uflag ? weno_size : 0) + offset;
      ww3 = w3 + (upw < 0 ? 2*weno_size : 0) + (uflag ? weno_size : 0) + offset;

      _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
      _ArrayIndexnD_(ndims,p,bounds_inter,indexC,0);
      _ArrayCopy1D_(indexC,indexI,ndims);

      if (upw > 0) {
          indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
          qm3 = qm1 - 2*stride;
          qm2 = qm1 -   stride;
          qp1 = qm1 +   stride;
          qp2 = qm1 + 2*stride;
      } else {
          indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
          qm3 = qm1 + 2*stride;
          qm2 = qm1 +   stride;
          qp1 = qm1 -   stride;
          qp2 = qm1 - 2*stride;
      }

      /* Defining stencil points */
      const double *fm3, *fm2, *fm1, *fp1, *fp2;
      fm3 = (fC+qm3*nvars);
      fm2 = (fC+qm2*nvars);
      fm1 = (fC+qm1*nvars);
      fp1 = (fC+qp1*nvars);
      fp2 = (fC+qp2*nvars);

      /* Candidate stencils and their optimal weights*/
      double f1[GPU_MAX_NVARS], f2[GPU_MAX_NVARS], f3[GPU_MAX_NVARS];

      _ArrayAXBYCZ_(f1,(2*GPU_ONE_SIXTH),fm3,(-7*GPU_ONE_SIXTH) ,fm2,(11*GPU_ONE_SIXTH) ,fm1,nvars);
      _ArrayAXBYCZ_(f2,(-GPU_ONE_SIXTH) ,fm2,(5*GPU_ONE_SIXTH)  ,fm1,(2*GPU_ONE_SIXTH)  ,fp1,nvars);
      _ArrayAXBYCZ_(f3,(2*GPU_ONE_SIXTH),fm1,(5*GPU_ONE_SIXTH)  ,fp1,(-GPU_ONE_SIXTH)   ,fp2,nvars);

      /* calculate WENO weights */
      const double *cur_w1, *cur_w2, *cur_w3;
      cur_w1 = (ww1+p*nvars);
      cur_w2 = (ww2+p*nvars);
      cur_w3 = (ww3+p*nvars);

      _ArrayMultiply3Add1D_((fI+p*nvars),cur_w1,f1,cur_w2,f2,cur_w3,f3,nvars);
    }

    return;
}

/*! @brief 5th order WENO reconstruction (component-wise) on a uniform grid

    Computes the interpolated values of the first primitive of a function \f${\bf f}\left({\bf u}\right)\f$
    at the interfaces from the cell-centered values of the function using the fifth order WENO scheme on a
    uniform grid. The first primitive is defined as a function \f${\bf h}\left({\bf u}\right)\f$ that satisfies:
    \f{equation}{
      {\bf f}\left({\bf u}\left(x\right)\right) = \frac{1}{\Delta x} \int_{x-\Delta x/2}^{x+\Delta x/2} {\bf h}\left({\bf u}\left(\zeta\right)\right)d\zeta,
    \f}
    where \f$x\f$ is the spatial coordinate along the dimension of the interpolation. This function computes the 5th order WENO numerical approximation
    \f$\hat{\bf f}_{j+1/2} \approx {\bf h}_{j+1/2}\f$ as the convex combination of three 3rd order methods:
    \f{align}{
        &\ \omega_1\ \times\ \left[ \hat{\bf f}_{j+1/2}^1 = \frac{1}{3} {\bf f}_{j-2} - \frac{7}{6} {\bf f}_{j-1} + \frac{11}{6} {\bf f}_j \right]\\
      + &\ \omega_2\ \times\ \left[ \hat{\bf f}_{j+1/2}^2 = -\frac{1}{6} {\bf f}_{j-1} + \frac{5}{6} {\bf f}_j + \frac{1}{3} {\bf f}_{j+1} \right]\\
      + &\ \omega_3\ \times\ \left[ \hat{\bf f}_{j+1/2}^3 = \frac{1}{3} {\bf f}_j + \frac{5}{6} {\bf f}_{j+1} - \frac{1}{6} {\bf f}_{j+2}  \right]\\
      \Rightarrow &\ \hat{\bf f}_{j+1/2} = \frac{\omega_1}{3} {\bf f}_{j-2} - \frac{1}{6}(7\omega_1+\omega_2){\bf f}_{j-1} + \frac{1}{6}(11\omega_1+5\omega_2+2\omega_3){\bf f}_j + \frac{1}{6}(2\omega_2+5\omega_3){\bf f}_{j+1} - \frac{\omega_3}{6}{\bf f}_{j+2},
    \f}
    where \f$\omega_k; k=1,2,3\f$ are the nonlinear WENO weights computed in WENOFifthOrderCalculateWeights() (note that the \f$\omega\f$ are different for each component of the vector \f$\hat{\bf f}\f$).

    \b Implementation \b Notes:
    + This method assumes a uniform grid in the spatial dimension corresponding to the interpolation.
    + The method described above corresponds to a left-biased interpolation. The corresponding right-biased
      interpolation can be obtained by reflecting the equations about interface j+1/2.
    + The scalar interpolation method is applied to the vector function in a component-wise manner.
    + The function computes the interpolant for the entire grid in one call. It loops over all the grid lines along the interpolation direction
      and carries out the 1D interpolation along these grid lines.
    + Location of cell-centers and cell interfaces along the spatial dimension of the interpolation is shown in the following figure:
      @image html chap1_1Ddomain.png
      @image latex chap1_1Ddomain.eps width=0.9\textwidth

    \b Function \b arguments:

    Argument  | Type      | Explanation
    --------- | --------- | ---------------------------------------------
    fI        | double*   | Array to hold the computed interpolant at the grid interfaces. This array must have the same layout as the solution, but with \b no \b ghost \b points. Its size should be the same as u in all dimensions, except dir (the dimension along which to interpolate) along which it should be larger by 1 (number of interfaces is 1 more than the number of interior cell centers).
    fC        | double*   | Array with the cell-centered values of the flux function \f${\bf f}\left({\bf u}\right)\f$. This array must have the same layout and size as the solution, \b with \b ghost \b points.
    u         | double*   | The solution array \f${\bf u}\f$ (with ghost points). If the interpolation is characteristic based, this is needed to compute the eigendecomposition. For a multidimensional problem, the layout is as follows: u is a contiguous 1D array of size (nvars*dim[0]*dim[1]*...*dim[D-1]) corresponding to the multi-dimensional solution, with the following ordering - nvars, dim[0], dim[1], ..., dim[D-1], where nvars is the number of solution components (#HyPar::nvars), dim is the local size (#HyPar::dim_local), D is the number of spatial dimensions.
    x         | double*   | The grid array (with ghost points). This is used only by non-uniform-grid interpolation methods. For multidimensional problems, the layout is as follows: x is a contiguous 1D array of size (dim[0]+dim[1]+...+dim[D-1]), with the spatial coordinates along dim[0] stored from 0,...,dim[0]-1, the spatial coordinates along dim[1] stored along dim[0],...,dim[0]+dim[1]-1, and so forth.
    upw       | int       | Upwinding direction: if positive, a left-biased interpolant will be computed; if negative, a right-biased interpolant will be computed. If the interpolation method is central, then this has no effect.
    dir       | int       | Spatial dimension along which to interpolate (eg: 0 for 1D; 0 or 1 for 2D; 0,1 or 2 for 3D)
    s         | void*     | Solver object of type #HyPar: the following variables are needed - #HyPar::ghosts, #HyPar::ndims, #HyPar::nvars, #HyPar::dim_local.
    m         | void*     | MPI object of type #MPIVariables: this is needed only by compact interpolation method that need to solve a global implicit system across MPI ranks.
    uflag     | int       | A flag indicating if the function being interpolated \f${\bf f}\f$ is the solution itself \f${\bf u}\f$ (if 1, \f${\bf f}\left({\bf u}\right) \equiv {\bf u}\f$).


    \b Reference:
    + Jiang, G.-S., Shu, C.-W., Efficient Implementation of Weighted ENO Schemes, J. Comput. Phys., 126 (1), 1996, pp. 202-228, http://dx.doi.org/10.1006/jcph.1996.0130
 */
extern "C" int gpuInterp1PrimFifthOrderWENO(
    double       *fI,  /*!< Array of interpolated function values at the interfaces */
    double *fC,  /*!< Array of cell-centered values of the function \f${\bf f}\left({\bf u}\right)\f$ */
    double *u,   /*!< Array of cell-centered values of the solution \f${\bf u}\f$ */
    double *x,   /*!< Grid coordinates */
    int    upw,  /*!< Upwind direction (left or right biased) */
    int    dir,  /*!< Spatial dimension along which to interpolation */
    void   *s,   /*!< Object of type #HyPar containing solver-related variables */
    void   *m,   /*!< Object of type #MPIVariables containing MPI-related variables */
    int    uflag /*!< Flag to indicate if \f$f(u) \equiv u\f$, i.e, if the solution is being reconstructed */
    )
{
    HyPar           *solver = (HyPar*)          s;
    WENOParameters  *weno   = (WENOParameters*) solver->interp;

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
    int nblocks = (npoints_grid - 1) / GPU_THREADS_PER_BLOCK + 1;

#if defined(GPU_STAT)
    cudaEvent_t start, stop;
    float milliseconds = 0;
    checkCuda( cudaEventCreate(&start) );
    checkCuda( cudaEventCreate(&stop) );

    int weno_memory_accessed = 3*npoints_grid*nvars*sizeof(double);
    int fI_memory_accessed = npoints_grid*nvars*sizeof(double);
    int fC_memory_accessed = 1;
    for (int d=0; d<ndims; d++) {
      if (d == dir) fC_memory_accessed *= (dim[d]+2*ghosts);
      else          fC_memory_accessed *= dim[d];
    }
    fC_memory_accessed *= nvars*sizeof(double);

    checkCuda( cudaEventRecord(start) );
#endif

    Interp1PrimFifthOrderWENO_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        npoints_grid, ndims, dir, ghosts, nvars, weno->size, offset, stride[dir], upw, uflag,
        solver->gpu_dim_local, fC, weno->w1, weno->w2, weno->w3, fI
    );
    cudaDeviceSynchronize();

#if defined(GPU_STAT)
    checkCuda( cudaEventRecord(stop) );
    checkCuda( cudaEventSynchronize(stop) );
    checkCuda( cudaEventElapsedTime(&milliseconds, start, stop) );

    printf("%-50s GPU time (secs) = %.6f dir = %d bandwidth (GB/s) = %6.2f\n",
            "Interp1PrimFifthOrderWENO", milliseconds*1e-3, dir,
            (1e-6*(weno_memory_accessed+fI_memory_accessed+fC_memory_accessed))/milliseconds);

    checkCuda( cudaEventDestroy(start) );
    checkCuda( cudaEventDestroy(stop) );
#endif

    return 0;
}

#else

/*! Kernel for gpuInterp1PrimFifthOrderWENO() */
__global__
void Interp1PrimFifthOrderWENO_kernel(
    int npoints_grid,
    int npoints_local_wghosts,
    int ndims,
    int dir,
    int ghosts,
    int nvars,
    int weno_size,
    int offset,
    int stride,
    int upw,
    int uflag,
    const int *dim,
    const double *fC,
    const double *w1,
    const double *w2,
    const double *w3,
    double *fI
)
{
    int p = threadIdx.x + (blockDim.x * blockIdx.x);
    if (p < npoints_grid) {
      int    bounds_inter[GPU_MAX_NDIMS], indexC[GPU_MAX_NDIMS], indexI[GPU_MAX_NDIMS];
      int qm1,qm2,qm3,qp1,qp2;
      const double *ww1, *ww2, *ww3;

      ww1 = w1 + (upw < 0 ? 2*weno_size : 0) + (uflag ? weno_size : 0) + offset;
      ww2 = w2 + (upw < 0 ? 2*weno_size : 0) + (uflag ? weno_size : 0) + offset;
      ww3 = w3 + (upw < 0 ? 2*weno_size : 0) + (uflag ? weno_size : 0) + offset;

      _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
      _ArrayIndexnD_(ndims,p,bounds_inter,indexC,0);
      _ArrayCopy1D_(indexC,indexI,ndims);

      if (upw > 0) {
          indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
          qm3 = qm1 - 2*stride;
          qm2 = qm1 -   stride;
          qp1 = qm1 +   stride;
          qp2 = qm1 + 2*stride;
      } else {
          indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
          qm3 = qm1 + 2*stride;
          qm2 = qm1 +   stride;
          qp1 = qm1 -   stride;
          qp2 = qm1 - 2*stride;
      }

      /* Defining stencil points */
      const double *fm3, *fm2, *fm1, *fp1, *fp2;
      /* Candidate stencils and their optimal weights*/
      double f1, f2, f3;

      int l = p;
      for (int j = 0; j < nvars; j++) {
        fm3 = (fC+qm3);
        fm2 = (fC+qm2);
        fm1 = (fC+qm1);
        fp1 = (fC+qp1);
        fp2 = (fC+qp2);

        f1 = (2*GPU_ONE_SIXTH)*fm3[0] + (-7*GPU_ONE_SIXTH)*fm2[0] + (11*GPU_ONE_SIXTH)*fm1[0];
        f2 = ( -GPU_ONE_SIXTH)*fm2[0] + ( 5*GPU_ONE_SIXTH)*fm1[0] + ( 2*GPU_ONE_SIXTH)*fp1[0];
        f3 = (2*GPU_ONE_SIXTH)*fm1[0] + ( 5*GPU_ONE_SIXTH)*fp1[0] + (  -GPU_ONE_SIXTH)*fp2[0];

        qm3 += npoints_local_wghosts;
        qm2 += npoints_local_wghosts;
        qm1 += npoints_local_wghosts;
        qp1 += npoints_local_wghosts;
        qp2 += npoints_local_wghosts;

        /*
        _ArrayAXBYCZ_(f1,(2*GPU_ONE_SIXTH),fm3,(-7*GPU_ONE_SIXTH) ,fm2,(11*GPU_ONE_SIXTH) ,fm1,nvars);
        _ArrayAXBYCZ_(f2,(-GPU_ONE_SIXTH) ,fm2,(5*GPU_ONE_SIXTH)  ,fm1,(2*GPU_ONE_SIXTH)  ,fp1,nvars);
        _ArrayAXBYCZ_(f3,(2*GPU_ONE_SIXTH),fm1,(5*GPU_ONE_SIXTH)  ,fp1,(-GPU_ONE_SIXTH)   ,fp2,nvars);
        */

        /* calculate WENO weights */
        fI[l] = (ww1+l)[0]*f1 + (ww2+l)[0]*f2 + (ww3+l)[0]*f3;
        l += npoints_grid;
        //_ArrayMultiply3Add1D_((fI+p),cur_w1,f1,cur_w2,f2,cur_w3,f3,1);
      }
    }

    return;
}

/*! @brief 5th order WENO reconstruction (component-wise) on a uniform grid

    Computes the interpolated values of the first primitive of a function \f${\bf f}\left({\bf u}\right)\f$
    at the interfaces from the cell-centered values of the function using the fifth order WENO scheme on a
    uniform grid. The first primitive is defined as a function \f${\bf h}\left({\bf u}\right)\f$ that satisfies:
    \f{equation}{
      {\bf f}\left({\bf u}\left(x\right)\right) = \frac{1}{\Delta x} \int_{x-\Delta x/2}^{x+\Delta x/2} {\bf h}\left({\bf u}\left(\zeta\right)\right)d\zeta,
    \f}
    where \f$x\f$ is the spatial coordinate along the dimension of the interpolation. This function computes the 5th order WENO numerical approximation
    \f$\hat{\bf f}_{j+1/2} \approx {\bf h}_{j+1/2}\f$ as the convex combination of three 3rd order methods:
    \f{align}{
        &\ \omega_1\ \times\ \left[ \hat{\bf f}_{j+1/2}^1 = \frac{1}{3} {\bf f}_{j-2} - \frac{7}{6} {\bf f}_{j-1} + \frac{11}{6} {\bf f}_j \right]\\
      + &\ \omega_2\ \times\ \left[ \hat{\bf f}_{j+1/2}^2 = -\frac{1}{6} {\bf f}_{j-1} + \frac{5}{6} {\bf f}_j + \frac{1}{3} {\bf f}_{j+1} \right]\\
      + &\ \omega_3\ \times\ \left[ \hat{\bf f}_{j+1/2}^3 = \frac{1}{3} {\bf f}_j + \frac{5}{6} {\bf f}_{j+1} - \frac{1}{6} {\bf f}_{j+2}  \right]\\
      \Rightarrow &\ \hat{\bf f}_{j+1/2} = \frac{\omega_1}{3} {\bf f}_{j-2} - \frac{1}{6}(7\omega_1+\omega_2){\bf f}_{j-1} + \frac{1}{6}(11\omega_1+5\omega_2+2\omega_3){\bf f}_j + \frac{1}{6}(2\omega_2+5\omega_3){\bf f}_{j+1} - \frac{\omega_3}{6}{\bf f}_{j+2},
    \f}
    where \f$\omega_k; k=1,2,3\f$ are the nonlinear WENO weights computed in WENOFifthOrderCalculateWeights() (note that the \f$\omega\f$ are different for each component of the vector \f$\hat{\bf f}\f$).

    \b Implementation \b Notes:
    + This method assumes a uniform grid in the spatial dimension corresponding to the interpolation.
    + The method described above corresponds to a left-biased interpolation. The corresponding right-biased
      interpolation can be obtained by reflecting the equations about interface j+1/2.
    + The scalar interpolation method is applied to the vector function in a component-wise manner.
    + The function computes the interpolant for the entire grid in one call. It loops over all the grid lines along the interpolation direction
      and carries out the 1D interpolation along these grid lines.
    + Location of cell-centers and cell interfaces along the spatial dimension of the interpolation is shown in the following figure:
      @image html chap1_1Ddomain.png
      @image latex chap1_1Ddomain.eps width=0.9\textwidth

    \b Function \b arguments:

    Argument  | Type      | Explanation
    --------- | --------- | ---------------------------------------------
    fI        | double*   | Array to hold the computed interpolant at the grid interfaces. This array must have the same layout as the solution, but with \b no \b ghost \b points. Its size should be the same as u in all dimensions, except dir (the dimension along which to interpolate) along which it should be larger by 1 (number of interfaces is 1 more than the number of interior cell centers).
    fC        | double*   | Array with the cell-centered values of the flux function \f${\bf f}\left({\bf u}\right)\f$. This array must have the same layout and size as the solution, \b with \b ghost \b points.
    u         | double*   | The solution array \f${\bf u}\f$ (with ghost points). If the interpolation is characteristic based, this is needed to compute the eigendecomposition. For a multidimensional problem, the layout is as follows: u is a contiguous 1D array of size (nvars*dim[0]*dim[1]*...*dim[D-1]) corresponding to the multi-dimensional solution, with the following ordering - nvars, dim[0], dim[1], ..., dim[D-1], where nvars is the number of solution components (#HyPar::nvars), dim is the local size (#HyPar::dim_local), D is the number of spatial dimensions.
    x         | double*   | The grid array (with ghost points). This is used only by non-uniform-grid interpolation methods. For multidimensional problems, the layout is as follows: x is a contiguous 1D array of size (dim[0]+dim[1]+...+dim[D-1]), with the spatial coordinates along dim[0] stored from 0,...,dim[0]-1, the spatial coordinates along dim[1] stored along dim[0],...,dim[0]+dim[1]-1, and so forth.
    upw       | int       | Upwinding direction: if positive, a left-biased interpolant will be computed; if negative, a right-biased interpolant will be computed. If the interpolation method is central, then this has no effect.
    dir       | int       | Spatial dimension along which to interpolate (eg: 0 for 1D; 0 or 1 for 2D; 0,1 or 2 for 3D)
    s         | void*     | Solver object of type #HyPar: the following variables are needed - #HyPar::ghosts, #HyPar::ndims, #HyPar::nvars, #HyPar::dim_local.
    m         | void*     | MPI object of type #MPIVariables: this is needed only by compact interpolation method that need to solve a global implicit system across MPI ranks.
    uflag     | int       | A flag indicating if the function being interpolated \f${\bf f}\f$ is the solution itself \f${\bf u}\f$ (if 1, \f${\bf f}\left({\bf u}\right) \equiv {\bf u}\f$).


    \b Reference:
    + Jiang, G.-S., Shu, C.-W., Efficient Implementation of Weighted ENO Schemes, J. Comput. Phys., 126 (1), 1996, pp. 202-228, http://dx.doi.org/10.1006/jcph.1996.0130
 */
extern "C" int gpuInterp1PrimFifthOrderWENO(
    double       *fI,  /*!< Array of interpolated function values at the interfaces */
    double *fC,  /*!< Array of cell-centered values of the function \f${\bf f}\left({\bf u}\right)\f$ */
    double *u,   /*!< Array of cell-centered values of the solution \f${\bf u}\f$ */
    double *x,   /*!< Grid coordinates */
    int    upw,  /*!< Upwind direction (left or right biased) */
    int    dir,  /*!< Spatial dimension along which to interpolation */
    void   *s,   /*!< Object of type #HyPar containing solver-related variables */
    void   *m,   /*!< Object of type #MPIVariables containing MPI-related variables */
    int    uflag /*!< Flag to indicate if \f$f(u) \equiv u\f$, i.e, if the solution is being reconstructed */
    )
{
    HyPar           *solver = (HyPar*)          s;
    WENOParameters  *weno   = (WENOParameters*) solver->interp;

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
    int nblocks = (npoints_grid - 1) / GPU_THREADS_PER_BLOCK + 1;

#if defined(GPU_STAT)
    cudaEvent_t start, stop;
    float milliseconds = 0;
    checkCuda( cudaEventCreate(&start) );
    checkCuda( cudaEventCreate(&stop) );


    int weno_memory_accessed = 3*npoints_grid*nvars*sizeof(double);
    int fI_memory_accessed = npoints_grid*nvars*sizeof(double);
    int fC_memory_accessed = 1;
    for (int d=0; d<ndims; d++) {
      if (d == dir) fC_memory_accessed *= (dim[d]+2*ghosts);
      else          fC_memory_accessed *= dim[d];
    }
    fC_memory_accessed *= nvars*sizeof(double);

    checkCuda( cudaEventRecord(start) );
#endif

    Interp1PrimFifthOrderWENO_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        npoints_grid, solver->npoints_local_wghosts, ndims, dir, ghosts, nvars, weno->size, offset, stride[dir], upw, uflag,
        solver->gpu_dim_local, fC, weno->w1, weno->w2, weno->w3, fI
    );
#if defined(GPU_STAT)
    checkCuda( cudaEventRecord(stop) );
    checkCuda( cudaEventSynchronize(stop) );
#endif

    cudaDeviceSynchronize();

#if defined(GPU_STAT)
    checkCuda( cudaEventElapsedTime(&milliseconds, start, stop) );

    printf("%-50s GPU time (secs) = %.6f dir = %d bandwidth (GB/s) = %6.2f\n",
            "Interp1PrimFifthOrderWENO2", milliseconds*1e-3, dir,
            (1e-6*(weno_memory_accessed+fI_memory_accessed+fC_memory_accessed))/milliseconds);

    checkCuda( cudaEventDestroy(start) );
    checkCuda( cudaEventDestroy(stop) );
#endif

    return 0;
}

#endif
