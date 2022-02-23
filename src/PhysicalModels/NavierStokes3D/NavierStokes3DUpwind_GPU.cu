/*! @file NavierStokes3DUpwind_GPU.cu
    @author Youngdae Kim
    @brief Contains functions to compute the upwind flux at grid interfaces for the 3D Navier Stokes equations.
*/
#include <math.h>
#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>
#include <mathfunctions.h>
#include <matmult_native.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

static const int dummy = 1;

#ifdef CUDA_VAR_ORDERDING_AOS

/* kernel for gpuNavierStokes3DUpwindRusanov() */
__global__
void gpuNavierStokes3DUpwindRusanov_kernel(
  int    npoints_grid,
  int    dir,
  int    ghosts,
  double gamma,
  const int    * __restrict__ dim,
  const double * __restrict__ grav_field_g,
  const double * __restrict__ fL,
  const double * __restrict__ fR,
  const double * __restrict__ uL,
  const double * __restrict__ uR,
  const double * __restrict__ u,
  double       * __restrict__ fI
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    double udiff[_MODEL_NVARS_],uavg[_MODEL_NVARS_];
    int bounds_inter[_MODEL_NDIMS_], index_inter[_MODEL_NDIMS_],
        indexL[_MODEL_NDIMS_], indexR[_MODEL_NDIMS_];

    bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[2] = dim[2]; bounds_inter[dir] += 1;
    _ArrayIndexnD_(_MODEL_NDIMS_,p,bounds_inter,index_inter,0);
    _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
    _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
    int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,ghosts,pL);
    int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,ghosts,pR);
    int q = p*_MODEL_NVARS_;

    /* Modified Rusanov's upwinding scheme */

    udiff[0] = 0.5 * (uR[q+0] - uL[q+0]);
    udiff[1] = 0.5 * (uR[q+1] - uL[q+1]);
    udiff[2] = 0.5 * (uR[q+2] - uL[q+2]);
    udiff[3] = 0.5 * (uR[q+3] - uL[q+3]);
    udiff[4] = 0.5 * (uR[q+4] - uL[q+4]);

    _NavierStokes3DRoeAverage_(uavg,_NavierStokes3D_stride_,(u+_MODEL_NVARS_*pL),(u+_MODEL_NVARS_*pR),gamma);

    double c, vel[_MODEL_NDIMS_], rho,E,P;
    _NavierStokes3DGetFlowVar_((u+_MODEL_NVARS_*pL),_NavierStokes3D_stride_,rho,vel[0],vel[1],vel[2],E,P,gamma);
    c = sqrt(gamma*P/rho);
    double alphaL = c + absolute(vel[dir]);
    _NavierStokes3DGetFlowVar_((u+_MODEL_NVARS_*pR),_NavierStokes3D_stride_,rho,vel[0],vel[1],vel[2],E,P,gamma);
    c = sqrt(gamma*P/rho);
    double alphaR = c + absolute(vel[dir]);
    _NavierStokes3DGetFlowVar_(uavg,dummy,rho,vel[0],vel[1],vel[2],E,P,gamma);
    c = sqrt(gamma*P/rho);
    double alphaavg = c + absolute(vel[dir]);

    double kappa  = max(grav_field_g[pL],grav_field_g[pR]);
    double alpha  = kappa*max3(alphaL,alphaR,alphaavg);

    fI[q+0] = 0.5*(fL[q+0]+fR[q+0])-alpha*udiff[0];
    fI[q+1] = 0.5*(fL[q+1]+fR[q+1])-alpha*udiff[1];
    fI[q+2] = 0.5*(fL[q+2]+fR[q+2])-alpha*udiff[2];
    fI[q+3] = 0.5*(fL[q+3]+fR[q+3])-alpha*udiff[3];
    fI[q+4] = 0.5*(fL[q+4]+fR[q+4])-alpha*udiff[4];
  }
  return;
}

/*! Rusanov's upwinding scheme.
    \f{equation}{
      {\bf f}_{j+1/2} = \frac{1}{2}\left[ {\bf f}_{j+1/2}^L + {\bf f}_{j+1/2}^R
                         - \max_{j,j+1} \nu_j \left( {\bf u}_{j+1/2}^R - {\bf u}_{j+1/2}^L  \right)\right]
    \f}
    where \f$\nu = c + \left|u\right|\f$.
    + Rusanov, V. V., "The calculation of the interaction of non-stationary shock waves and obstacles," USSR
    Computational Mathematics and Mathematical Physics, Vol. 1, No. 2, 1962, pp. 304–320

    This upwinding scheme is modified for the balanced discretization of the 3D Navier Stokes equations when
    there is a non-zero gravitational force. See the reference below. For flows without any gravitational forces,
    it reduces to its original form.
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580

*/
extern "C" int gpuNavierStokes3DUpwindRusanov(
  double  *fI, /*!< Computed upwind interface flux */
  double  *fL, /*!< Left-biased reconstructed interface flux */
  double  *fR, /*!< Right-biased reconstructed interface flux */
  double  *uL, /*!< Left-biased reconstructed interface solution */
  double  *uR, /*!< Right-biased reconstructed interface solution */
  double  *u,  /*!< Cell-centered solution */
  int     dir, /*!< Spatial dimension (x,y, or z) */
  void    *s,  /*!< Solver object of type #HyPar */
  double  t    /*!< Current solution time */
)
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;
  int             *dim    = solver->dim_local;

  int bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;
  int npoints_grid; _ArrayProduct1D_(bounds_inter,_MODEL_NDIMS_,npoints_grid);
  int nblocks = (npoints_grid-1)/GPU_THREADS_PER_BLOCK + 1;

  gpuNavierStokes3DUpwindRusanov_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    npoints_grid, dir, solver->ghosts, param->gamma, solver->gpu_dim_local,
    param->gpu_grav_field_g, fL, fR, uL, uR, u, fI
  );
  cudaDeviceSynchronize();

  return 0;
}

/* kernel for gpuNavierStokes3DUpwindRoe() */
__global__
void gpuNavierStokes3DUpwindRoe_kernel(
  int    npoints_grid,
  int    dir,
  int    ghosts,
  double gamma,
  const int    * __restrict__ dim,
  const double * __restrict__ grav_field_g,
  const double * __restrict__ fL,
  const double * __restrict__ fR,
  const double * __restrict__ uL,
  const double * __restrict__ uR,
  const double * __restrict__ u,
  double       * __restrict__ fI
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_],
           L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_],
           modA[_MODEL_NVARS_*_MODEL_NVARS_];
    double udiff[_MODEL_NVARS_],uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];
    int bounds_inter[_MODEL_NDIMS_], index_inter[_MODEL_NDIMS_],
        indexL[_MODEL_NDIMS_], indexR[_MODEL_NDIMS_];

    bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[2] = dim[2]; bounds_inter[dir] += 1;
    _ArrayIndexnD_(_MODEL_NDIMS_,p,bounds_inter,index_inter,0);
    _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
    _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
    int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,ghosts,pL);
    int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,ghosts,pR);

    p *= _MODEL_NVARS_;

    /* Modified Roe's upwinding scheme */

    udiff[0] = 0.5 * (uR[p+0] - uL[p+0]);
    udiff[1] = 0.5 * (uR[p+1] - uL[p+1]);
    udiff[2] = 0.5 * (uR[p+2] - uL[p+2]);
    udiff[3] = 0.5 * (uR[p+3] - uL[p+3]);
    udiff[4] = 0.5 * (uR[p+4] - uL[p+4]);

    _NavierStokes3DRoeAverage_        (uavg,_NavierStokes3D_stride_,(u+_MODEL_NVARS_*pL),(u+_MODEL_NVARS_*pR),gamma);
    _NavierStokes3DEigenvalues_       (uavg,dummy,D,gamma,dir);
    _NavierStokes3DLeftEigenvectors_  (uavg,dummy,L,gamma,dir);
    _NavierStokes3DRightEigenvectors_ (uavg,dummy,R,gamma,dir);

    /* Harten's Entropy Fix - Page 362 of Leveque */
    int k;
    double delta = 0.000001, delta2 = delta*delta;
    k=0;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
    k=6;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
    k=12; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
    k=18; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
    k=24; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );

    MatMult5(_MODEL_NVARS_,DL,D,L);
    MatMult5(_MODEL_NVARS_,modA,R,DL);
    MatVecMult5(_MODEL_NVARS_,udiss,modA,udiff);

    fI[p+0] = 0.5 * (fL[p+0]+fR[p+0]) - udiss[0];
    fI[p+1] = 0.5 * (fL[p+1]+fR[p+1]) - udiss[1];
    fI[p+2] = 0.5 * (fL[p+2]+fR[p+2]) - udiss[2];
    fI[p+3] = 0.5 * (fL[p+3]+fR[p+3]) - udiss[3];
    fI[p+4] = 0.5 * (fL[p+4]+fR[p+4]) - udiss[4];
  }
  return;
}

/*! Roe's upwinding scheme.
    \f{equation}{
      {\bf f}_{j+1/2} = \frac{1}{2}\left[ {\bf f}_{j+1/2}^L + {\bf f}_{j+1/2}^R
                         - \left| A\left({\bf u}_{j+1/2}^L,{\bf u}_{j+1/2}^R\right) \right|
                           \left( {\bf u}_{j+1/2}^R - {\bf u}_{j+1/2}^L  \right)\right]
    \f}
    + Roe, P. L., “Approximate Riemann solvers, parameter vectors, and difference schemes,” Journal of
    Computational Physics, Vol. 43, No. 2, 1981, pp. 357–372, http://dx.doi.org/10.1016/0021-9991(81)90128-5.

    This upwinding scheme is modified for the balanced discretization of the 3D Navier Stokes equations when
    there is a non-zero gravitational force. See the reference below. For flows without any gravitational forces,
    it reduces to its original form.
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580

*/
extern "C" int gpuNavierStokes3DUpwindRoe(
  double  * __restrict__ fI, /*!< Computed upwind interface flux */
  double  * __restrict__ fL, /*!< Left-biased reconstructed interface flux */
  double  * __restrict__ fR, /*!< Right-biased reconstructed interface flux */
  double  * __restrict__ uL, /*!< Left-biased reconstructed interface solution */
  double  * __restrict__ uR, /*!< Right-biased reconstructed interface solution */
  double  * __restrict__ u,  /*!< Cell-centered solution */
  int     dir, /*!< Spatial dimension (x,y, or z) */
  void    *s,  /*!< Solver object of type #HyPar */
  double  t    /*!< Current solution time */
)
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;
  int             *dim    = solver->dim_local;

  int bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;
  int npoints_grid; _ArrayProduct1D_(bounds_inter,_MODEL_NDIMS_,npoints_grid);
  int nblocks = (npoints_grid-1)/GPU_THREADS_PER_BLOCK + 1;

  gpuNavierStokes3DUpwindRoe_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    npoints_grid, dir, solver->ghosts, param->gamma, solver->gpu_dim_local,
    param->gpu_grav_field_g, fL, fR, uL, uR, u, fI
  );
  cudaDeviceSynchronize();

  return 0;
}

#else

/* kernel for gpuNavierStokes3DUpwindRusanov() */
__global__
void gpuNavierStokes3DUpwindRusanov_kernel(
  int    npoints_grid,
  int    npoints_local_wghosts,
  int    dir,
  int    ghosts,
  double gamma,
  const int    * __restrict__ dim,
  const double * __restrict__ grav_field_g,
  const double * __restrict__ fL,
  const double * __restrict__ fR,
  const double * __restrict__ uL,
  const double * __restrict__ uR,
  const double * __restrict__ u,
  double       * __restrict__ fI
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    double udiff[_MODEL_NVARS_],uavg[_MODEL_NVARS_];
    int bounds_inter[_MODEL_NDIMS_], index_inter[_MODEL_NDIMS_],
        indexL[_MODEL_NDIMS_], indexR[_MODEL_NDIMS_];

    bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[2] = dim[2]; bounds_inter[dir] += 1;
    _ArrayIndexnD_(_MODEL_NDIMS_,p,bounds_inter,index_inter,0);
    _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
    _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
    int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,ghosts,pL);
    int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,ghosts,pR);

    /* Modified Rusanov's upwinding scheme */

    udiff[0] = 0.5 * (uR[p+             0] - uL[p+             0]);
    udiff[1] = 0.5 * (uR[p+  npoints_grid] - uL[p+  npoints_grid]);
    udiff[2] = 0.5 * (uR[p+2*npoints_grid] - uL[p+2*npoints_grid]);
    udiff[3] = 0.5 * (uR[p+3*npoints_grid] - uL[p+3*npoints_grid]);
    udiff[4] = 0.5 * (uR[p+4*npoints_grid] - uL[p+4*npoints_grid]);

    _NavierStokes3DRoeAverage_(uavg,npoints_local_wghosts,(u+pL),(u+pR),gamma);

    double c, vel[_MODEL_NDIMS_], rho,E,P;
    _NavierStokes3DGetFlowVar_((u+pL),npoints_local_wghosts,rho,vel[0],vel[1],vel[2],E,P,gamma);
    c = sqrt(gamma*P/rho);
    double alphaL = c + absolute(vel[dir]);
    _NavierStokes3DGetFlowVar_((u+pR),npoints_local_wghosts,rho,vel[0],vel[1],vel[2],E,P,gamma);
    c = sqrt(gamma*P/rho);
    double alphaR = c + absolute(vel[dir]);
    _NavierStokes3DGetFlowVar_(uavg,dummy,rho,vel[0],vel[1],vel[2],E,P,gamma);
    c = sqrt(gamma*P/rho);
    double alphaavg = c + absolute(vel[dir]);

    double kappa  = max(grav_field_g[pL],grav_field_g[pR]);
    double alpha  = kappa*max3(alphaL,alphaR,alphaavg);

    fI[p+0]              = 0.5*(fL[p+             0]+fR[p+             0])-alpha*udiff[0];
    fI[p+  npoints_grid] = 0.5*(fL[p+  npoints_grid]+fR[p+  npoints_grid])-alpha*udiff[1];
    fI[p+2*npoints_grid] = 0.5*(fL[p+2*npoints_grid]+fR[p+2*npoints_grid])-alpha*udiff[2];
    fI[p+3*npoints_grid] = 0.5*(fL[p+3*npoints_grid]+fR[p+3*npoints_grid])-alpha*udiff[3];
    fI[p+4*npoints_grid] = 0.5*(fL[p+4*npoints_grid]+fR[p+4*npoints_grid])-alpha*udiff[4];
  }
  return;
}

/*! Rusanov's upwinding scheme.
    \f{equation}{
      {\bf f}_{j+1/2} = \frac{1}{2}\left[ {\bf f}_{j+1/2}^L + {\bf f}_{j+1/2}^R
                         - \max_{j,j+1} \nu_j \left( {\bf u}_{j+1/2}^R - {\bf u}_{j+1/2}^L  \right)\right]
    \f}
    where \f$\nu = c + \left|u\right|\f$.
    + Rusanov, V. V., "The calculation of the interaction of non-stationary shock waves and obstacles," USSR
    Computational Mathematics and Mathematical Physics, Vol. 1, No. 2, 1962, pp. 304–320

    This upwinding scheme is modified for the balanced discretization of the 3D Navier Stokes equations when
    there is a non-zero gravitational force. See the reference below. For flows without any gravitational forces,
    it reduces to its original form.
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580

*/
extern "C" int gpuNavierStokes3DUpwindRusanov(
  double  *fI, /*!< Computed upwind interface flux */
  double  *fL, /*!< Left-biased reconstructed interface flux */
  double  *fR, /*!< Right-biased reconstructed interface flux */
  double  *uL, /*!< Left-biased reconstructed interface solution */
  double  *uR, /*!< Right-biased reconstructed interface solution */
  double  *u,  /*!< Cell-centered solution */
  int     dir, /*!< Spatial dimension (x,y, or z) */
  void    *s,  /*!< Solver object of type #HyPar */
  double  t    /*!< Current solution time */
)
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;
  int             *dim    = solver->dim_local;

  int bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;
  int npoints_grid; _ArrayProduct1D_(bounds_inter,_MODEL_NDIMS_,npoints_grid);
  int nblocks = (npoints_grid-1)/GPU_THREADS_PER_BLOCK + 1;

  gpuNavierStokes3DUpwindRusanov_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    npoints_grid, solver->npoints_local_wghosts, dir, solver->ghosts, param->gamma, solver->gpu_dim_local,
    param->gpu_grav_field_g, fL, fR, uL, uR, u, fI
  );
  cudaDeviceSynchronize();

  return 0;
}

/* kernel for gpuNavierStokes3DUpwindRoe() */
__global__
void gpuNavierStokes3DUpwindRoe_kernel(
  int    npoints_grid,
  int    npoints_local_wghosts,
  int    dir,
  int    ghosts,
  double gamma,
  const int    * __restrict__ dim,
  const double * __restrict__ grav_field_g,
  const double * __restrict__ fL,
  const double * __restrict__ fR,
  const double * __restrict__ uL,
  const double * __restrict__ uR,
  const double * __restrict__ u,
  double       * __restrict__ fI
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_],
           L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_],
           modA[_MODEL_NVARS_*_MODEL_NVARS_];
    double udiff[_MODEL_NVARS_],uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];
    int bounds_inter[_MODEL_NDIMS_], index_inter[_MODEL_NDIMS_],
        indexL[_MODEL_NDIMS_], indexR[_MODEL_NDIMS_];

    bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[2] = dim[2]; bounds_inter[dir] += 1;
    _ArrayIndexnD_(_MODEL_NDIMS_,p,bounds_inter,index_inter,0);
    _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
    _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
    int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,ghosts,pL);
    int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,ghosts,pR);

    /* Modified Roe's upwinding scheme */

    udiff[0] = 0.5 * (uR[p+0] - uL[p+0]);
    udiff[1] = 0.5 * (uR[p+  npoints_grid] - uL[p+  npoints_grid]);
    udiff[2] = 0.5 * (uR[p+2*npoints_grid] - uL[p+2*npoints_grid]);
    udiff[3] = 0.5 * (uR[p+3*npoints_grid] - uL[p+3*npoints_grid]);
    udiff[4] = 0.5 * (uR[p+4*npoints_grid] - uL[p+4*npoints_grid]);

    _NavierStokes3DRoeAverage_        (uavg,npoints_local_wghosts,(u+pL),(u+pR),gamma);
    _NavierStokes3DEigenvalues_       (uavg,dummy,D,gamma,dir);
    _NavierStokes3DLeftEigenvectors_  (uavg,dummy,L,gamma,dir);
    _NavierStokes3DRightEigenvectors_ (uavg,dummy,R,gamma,dir);

    /* Harten's Entropy Fix - Page 362 of Leveque */
    int k;
    double delta = 0.000001, delta2 = delta*delta;
    k=0;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
    k=6;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
    k=12; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
    k=18; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
    k=24; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );

    MatMult5(_MODEL_NVARS_,DL,D,L);
    MatMult5(_MODEL_NVARS_,modA,R,DL);
    MatVecMult5(_MODEL_NVARS_,udiss,modA,udiff);

    fI[p+0]              = 0.5 * (fL[p+0]+fR[p+0]) - udiss[0];
    fI[p+  npoints_grid] = 0.5 * (fL[p+  npoints_grid]+fR[p+  npoints_grid]) - udiss[1];
    fI[p+2*npoints_grid] = 0.5 * (fL[p+2*npoints_grid]+fR[p+2*npoints_grid]) - udiss[2];
    fI[p+3*npoints_grid] = 0.5 * (fL[p+3*npoints_grid]+fR[p+3*npoints_grid]) - udiss[3];
    fI[p+4*npoints_grid] = 0.5 * (fL[p+4*npoints_grid]+fR[p+4*npoints_grid]) - udiss[4];
  }
  return;
}

/*! Roe's upwinding scheme.
    \f{equation}{
      {\bf f}_{j+1/2} = \frac{1}{2}\left[ {\bf f}_{j+1/2}^L + {\bf f}_{j+1/2}^R
                         - \left| A\left({\bf u}_{j+1/2}^L,{\bf u}_{j+1/2}^R\right) \right|
                           \left( {\bf u}_{j+1/2}^R - {\bf u}_{j+1/2}^L  \right)\right]
    \f}
    + Roe, P. L., “Approximate Riemann solvers, parameter vectors, and difference schemes,” Journal of
    Computational Physics, Vol. 43, No. 2, 1981, pp. 357–372, http://dx.doi.org/10.1016/0021-9991(81)90128-5.

    This upwinding scheme is modified for the balanced discretization of the 3D Navier Stokes equations when
    there is a non-zero gravitational force. See the reference below. For flows without any gravitational forces,
    it reduces to its original form.
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580

*/
extern "C" int gpuNavierStokes3DUpwindRoe(
  double  * __restrict__ fI, /*!< Computed upwind interface flux */
  double  * __restrict__ fL, /*!< Left-biased reconstructed interface flux */
  double  * __restrict__ fR, /*!< Right-biased reconstructed interface flux */
  double  * __restrict__ uL, /*!< Left-biased reconstructed interface solution */
  double  * __restrict__ uR, /*!< Right-biased reconstructed interface solution */
  double  * __restrict__ u,  /*!< Cell-centered solution */
  int     dir, /*!< Spatial dimension (x,y, or z) */
  void    *s,  /*!< Solver object of type #HyPar */
  double  t    /*!< Current solution time */
)
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;
  int             *dim    = solver->dim_local;

  int bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;
  int npoints_grid; _ArrayProduct1D_(bounds_inter,_MODEL_NDIMS_,npoints_grid);
  int nblocks = (npoints_grid-1)/GPU_THREADS_PER_BLOCK + 1;

#if defined(GPU_STAT)
  cudaEvent_t startEvent, stopEvent;
  float milliseconds = 0;

  checkCuda( cudaEventCreate(&startEvent) );
  checkCuda( cudaEventCreate(&stopEvent) );
  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  gpuNavierStokes3DUpwindRoe_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    npoints_grid, solver->npoints_local_wghosts, dir, solver->ghosts, param->gamma, solver->gpu_dim_local,
    param->gpu_grav_field_g, fL, fR, uL, uR, u, fI
  );

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
#endif

  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );
  printf("%-50s GPU time (secs) = %.6f dir = %d\n",
          "NavierStokes3DUpwindRoe2", milliseconds*1e-3, dir);
  checkCuda( cudaEventDestroy(startEvent) );
  checkCuda( cudaEventDestroy(stopEvent) );
#endif




  return 0;
}

#endif
