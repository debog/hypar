/*! @file NavierStokes3DModifiedSolution_GPU.cu
    @author Youngdae Kim
    @brief Compute the modified solution for the well-balanced treatment of gravitational source terms.
*/
#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>
#include <physicalmodels/navierstokes3d.h>
#include <mpivars.h>
#include <hypar.h>

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel for gpuNavierStokes3DModifiedSolution() */
__global__
void gpuNavierStokes3DModifiedSolution_kernel(
  int    npoints_grid,
  double gamma,
  double inv_gamma_m1,
  const double * __restrict__ grav_field_f,
  const double * __restrict__ grav_field_g,
  const double * __restrict__ u,
  double       * __restrict__ uC
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    double rho, uvel, vvel, wvel, E, P;

    _NavierStokes3DGetFlowVar_((u+_MODEL_NVARS_*p),_NavierStokes3D_stride_,rho,uvel,vvel,wvel,E,P,gamma);
    uC[_MODEL_NVARS_*p+0] = u[_MODEL_NVARS_*p+0] * grav_field_f[p];
    uC[_MODEL_NVARS_*p+1] = u[_MODEL_NVARS_*p+1] * grav_field_f[p];
    uC[_MODEL_NVARS_*p+2] = u[_MODEL_NVARS_*p+2] * grav_field_f[p];
    uC[_MODEL_NVARS_*p+3] = u[_MODEL_NVARS_*p+3] * grav_field_f[p];
    uC[_MODEL_NVARS_*p+4] = (P*inv_gamma_m1)*(1.0/grav_field_g[p]) + (0.5*rho*(uvel*uvel+vvel*vvel+wvel*wvel))*grav_field_f[p];
  }
  return;
}

/*!
    This function computes the modified solution for the well-balanced treatment of the
    gravitational source terms. The modified solution vector is given by
    \f{equation}{
      {\bf u}^* = \left[\begin{array}{c} \rho \varrho^{-1}\left(x,y\right) \\ \rho u \varrho^{-1}\left(x,y\right) \\ \rho v \varrho^{-1}\left(x,y\right) \\ \rho w \varrho^{-1}\left(x,y\right) \\ e^* \end{array}\right]
    \f}
    where
    \f{equation}{
      e^* = \frac{p \varphi^{-1}\left(x,y\right)}{\gamma-1} + \frac{1}{2}\rho \varrho^{-1}\left(x,y\right) \left(u^2+v^2+w^2\right)
    \f}
    \f$\varrho\f$ and \f$\varphi\f$ are computed in #NavierStokes3DGravityField(). For flows without gravity, \f${\bf u}^* = {\bf u}\f$.

    References:
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, AIAA Journal, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580.
*/
extern "C" int gpuNavierStokes3DModifiedSolution(
  double  * __restrict__ uC,  /*!< Array to hold the computed modified solution */
  double  * __restrict__ u,   /*!< Solution vector array */
  int     d,    /*!< spatial dimension (not used) */
  void    * __restrict__ s,   /*!< Solver object of type #HyPar */
  void    * __restrict__ m,   /*!< MPI object of time #MPIVariables */
  double waqt   /*!< Current simulation time */
)
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;

  double inv_gamma_m1 = 1.0 / (param->gamma-1.0);
  int nblocks = (solver->npoints_local_wghosts-1)/GPU_THREADS_PER_BLOCK + 1;

  gpuNavierStokes3DModifiedSolution_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    solver->npoints_local_wghosts, param->gamma, inv_gamma_m1,
    param->gpu_grav_field_f, param->gpu_grav_field_g, u, uC
  );
  cudaDeviceSynchronize();

  return 0;
}

#else

/*! Kernel for gpuNavierStokes3DModifiedSolution() */
__global__
void gpuNavierStokes3DModifiedSolution_kernel(
  int    npoints_grid,
  double gamma,
  double inv_gamma_m1,
  const double * __restrict__ grav_field_f,
  const double * __restrict__ grav_field_g,
  const double * __restrict__ u,
  double       * __restrict__ uC
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    double rho, uvel, vvel, wvel, E, P;

    _NavierStokes3DGetFlowVar_((u+p),npoints_grid,rho,uvel,vvel,wvel,E,P,gamma);
    uC[p+0] = u[p+0] * grav_field_f[p];
    uC[p+  npoints_grid] = u[p+  npoints_grid] * grav_field_f[p];
    uC[p+2*npoints_grid] = u[p+2*npoints_grid] * grav_field_f[p];
    uC[p+3*npoints_grid] = u[p+3*npoints_grid] * grav_field_f[p];
    uC[p+4*npoints_grid] = (P*inv_gamma_m1)*(1.0/grav_field_g[p]) + (0.5*rho*(uvel*uvel+vvel*vvel+wvel*wvel))*grav_field_f[p];
  }
  return;
}

/*!
    This function computes the modified solution for the well-balanced treatment of the
    gravitational source terms. The modified solution vector is given by
    \f{equation}{
      {\bf u}^* = \left[\begin{array}{c} \rho \varrho^{-1}\left(x,y\right) \\ \rho u \varrho^{-1}\left(x,y\right) \\ \rho v \varrho^{-1}\left(x,y\right) \\ \rho w \varrho^{-1}\left(x,y\right) \\ e^* \end{array}\right]
    \f}
    where
    \f{equation}{
      e^* = \frac{p \varphi^{-1}\left(x,y\right)}{\gamma-1} + \frac{1}{2}\rho \varrho^{-1}\left(x,y\right) \left(u^2+v^2+w^2\right)
    \f}
    \f$\varrho\f$ and \f$\varphi\f$ are computed in #NavierStokes3DGravityField(). For flows without gravity, \f${\bf u}^* = {\bf u}\f$.

    References:
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, AIAA Journal, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580.
*/
extern "C" int gpuNavierStokes3DModifiedSolution(
  double  * __restrict__ uC,  /*!< Array to hold the computed modified solution */
  double  * __restrict__ u,   /*!< Solution vector array */
  int     d,                  /*!< spatial dimension (not used) */
  void    * __restrict__ s,   /*!< Solver object of type #HyPar */
  void    * __restrict__ m,   /*!< MPI object of time #MPIVariables */
  double waqt                 /*!< Current simulation time */
)
{
  HyPar           *solver = (HyPar*)          s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;

  double inv_gamma_m1 = 1.0 / (param->gamma-1.0);
  int nblocks = (solver->npoints_local_wghosts-1)/GPU_THREADS_PER_BLOCK + 1;

  gpuNavierStokes3DModifiedSolution_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    solver->npoints_local_wghosts, param->gamma, inv_gamma_m1,
    param->gpu_grav_field_f, param->gpu_grav_field_g, u, uC
  );
  cudaDeviceSynchronize();

  return 0;
}
#endif

