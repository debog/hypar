/*! @file NavierStokes3DFlux_GPU.cu
    @author Youngdae Kim
    @brief Functions to compute the hyperbolic flux for 3D Navier-Stokes equations
*/

#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>
#include <mathfunctions.h>
#include <matmult_native.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel for gpuNavierStokes3DFlux() */
__global__
void gpuNavierStokes3DFlux_kernel(
  int npoints_grid,
  int dir,
  double gamma,
  const double * __restrict__ u,
  double       * __restrict__ f
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    double rho, vx, vy, vz, e, P;
    _NavierStokes3DGetFlowVar_((u+_MODEL_NVARS_*p),_NavierStokes3D_stride_,rho,vx,vy,vz,e,P,gamma);
    _NavierStokes3DSetFlux_((f+_MODEL_NVARS_*p),_NavierStokes3D_stride_,rho,vx,vy,vz,e,P,dir);
  }
  return;
}

/*!
  Compute the hyperbolic flux function for the 3D Navier-Stokes equations:
  \f{eqnarray}{
    dir = x, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho u \\ \rho u^2 + p \\ \rho u v \\ \rho u w \\ (e+p)u \end{array}\right], \\
    dir = y, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho v \\ \rho u v \\ \rho v^2 + p \\ \rho v w \\ (e+p)v \end{array}\right], \\
    dir = z, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho w \\ \rho u w \\ \rho v w \\ \rho w^2 + p \\ (e+p)w \end{array}\right]
  \f}
  Note: the flux function needs to be computed at the ghost points as well.
*/
extern "C" int gpuNavierStokes3DFlux(
  double *__restrict__ f,  /*!< Array to hold the computed flux vector (same layout as u) */
  double *__restrict__ u,  /*!< Array with the solution vector */
  int    dir, /*!< Spatial dimension (x, y, or z) for which to compute the flux */
  void   *__restrict__ s,  /*!< Solver object of type #HyPar */
  double t  /*!< Current simulation time */
)
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes3D    *param  = (NavierStokes3D*) solver->physics;

  int nblocks = (solver->npoints_local_wghosts-1)/GPU_THREADS_PER_BLOCK + 1;

  gpuNavierStokes3DFlux_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    solver->npoints_local_wghosts, dir, param->gamma, u, f
  );
  cudaDeviceSynchronize();

  return 0;
}

#else

/*! Kernel for gpuNavierStokes3DFlux() */
__global__
void gpuNavierStokes3DFlux_kernel(
  int npoints_grid,
  int dir,
  double gamma,
  const double * __restrict__ u,
  double       * __restrict__ f
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    double rho, vx, vy, vz, e, P;
    _NavierStokes3DGetFlowVar_((u+p),npoints_grid,rho,vx,vy,vz,e,P,gamma);
    _NavierStokes3DSetFlux_((f+p),npoints_grid,rho,vx,vy,vz,e,P,dir);
  }
  return;
}

/*!
  Compute the hyperbolic flux function for the 3D Navier-Stokes equations:
  \f{eqnarray}{
    dir = x, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho u \\ \rho u^2 + p \\ \rho u v \\ \rho u w \\ (e+p)u \end{array}\right], \\
    dir = y, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho v \\ \rho u v \\ \rho v^2 + p \\ \rho v w \\ (e+p)v \end{array}\right], \\
    dir = z, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho w \\ \rho u w \\ \rho v w \\ \rho w^2 + p \\ (e+p)w \end{array}\right]
  \f}
  Note: the flux function needs to be computed at the ghost points as well.
*/
extern "C" int gpuNavierStokes3DFlux(
  double *__restrict__ f,  /*!< Array to hold the computed flux vector (same layout as u) */
  double *__restrict__ u,  /*!< Array with the solution vector */
  int    dir, /*!< Spatial dimension (x, y, or z) for which to compute the flux */
  void   *__restrict__ s,  /*!< Solver object of type #HyPar */
  double t  /*!< Current simulation time */
)
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes3D    *param  = (NavierStokes3D*) solver->physics;

  int nblocks = (solver->npoints_local_wghosts-1)/GPU_THREADS_PER_BLOCK + 1;

  gpuNavierStokes3DFlux_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    solver->npoints_local_wghosts, dir, param->gamma, u, f
  );
  cudaDeviceSynchronize();

  return 0;
}

#endif
