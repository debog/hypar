/*! @file NavierStokes2DFlux_GPU.cu
    @author Youngdae Kim
    @brief Functions to compute the hyperbolic flux for 2D Navier-Stokes equations
*/

#include <stdlib.h>
#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>
#include <mathfunctions.h>
#include <matmult_native.h>
#include <physicalmodels/navierstokes2d.h>
#include <hypar.h>

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel for gpuNavierStokes2DFlux() */
__global__
void NavierStokes2DFlux_kernel(
  int ngrid_points,
  int dir,
  double gamma,
  const double *u,
  double *f
)
{
  int p = threadIdx.x + (blockDim.x * blockIdx.x);
  if (p < ngrid_points) {
    double rho, vx, vy, e, P;
    _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*p),rho,vx,vy,e,P,gamma);
    _NavierStokes2DSetFlux_((f+_MODEL_NVARS_*p),rho,vx,vy,e,P,dir);
  }
  return;
}

/*!
  Compute the hyperbolic flux function for the 2D Navier-Stokes equations:
  \f{eqnarray}{
    dir = x, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho u \\ \rho u^2 + p \\ \rho u v \\ (e+p)u \end{array}\right], \\
    dir = y, & {\bf f}\left({\bf u}\right) = \left[\begin{array}{c} \rho v \\ \rho u v \\ \rho v^2 + p \\ (e+p)v \end{array}\right]
  \f}
  Note: the flux function needs to be computed at the ghost points as well.
*/
extern "C" int gpuNavierStokes2DFlux(
  double  *f, /*!< Array to hold the computed flux vector (same layout as u) */
  double  *u, /*!< Array with the solution vector */
  int     dir,   /*!< Spatial dimension (x or y) for which to compute the flux */
  void    *s,    /*!< Solver object of type #HyPar */
  double  t      /*!< Current simulation time */
)
{
  HyPar           *solver = (HyPar*)   s;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;
  int             *dim    = solver->dim_local;
  int             ghosts  = solver->ghosts;

  int ngrid_points = 1; for (int i=0; i<_MODEL_NDIMS_; i++) ngrid_points *= (dim[i]+2*ghosts);
  int nblocks = (ngrid_points - 1) / GPU_THREADS_PER_BLOCK + 1;

  clock_t cpu_start, cpu_end;
  double cpu_time = 0.0;

  cpu_start = clock();
  NavierStokes2DFlux_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    ngrid_points,dir,param->gamma,u,f
  );
  cudaDeviceSynchronize();
  cpu_end = clock();
  cpu_time += (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;

  return(0);
}

#endif
