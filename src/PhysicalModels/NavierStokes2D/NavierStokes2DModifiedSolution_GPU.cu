/*! @file NavierStokes2DModifiedSolution_GPU.cu
    @author Youngdae Kim
    @brief Compute the modified solution for the well-balanced treatment of gravitational source terms.
*/
#include <basic_gpu.h>
#include <arrayfunctions.h>
#include <physicalmodels/navierstokes2d.h>
#include <mpivars.h>
#include <hypar.h>

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel for gpuNavierStokes2DModifiedSolution() */
__global__
void NavierStokes2DModifiedSolution_kernel(
    int size,
    double gamma,
    double inv_gamma_m1,
    const double *grav_field_f,
    const double *grav_field_g,
    const double *u,
    double *uC
)
{
    int p = threadIdx.x + (blockDim.x * blockIdx.x);

    if (p < size) {
        double rho, uvel, vvel, E, P;
        _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*p),rho,uvel,vvel,E,P,gamma);
        uC[_MODEL_NVARS_*p+0] = u[_MODEL_NVARS_*p+0] * grav_field_f[p];
        uC[_MODEL_NVARS_*p+1] = u[_MODEL_NVARS_*p+1] * grav_field_f[p];
        uC[_MODEL_NVARS_*p+2] = u[_MODEL_NVARS_*p+2] * grav_field_f[p];
        uC[_MODEL_NVARS_*p+3] = (P*inv_gamma_m1)*(1.0/grav_field_g[p]) + (0.5*rho*(uvel*uvel+vvel*vvel))*grav_field_f[p];
    }

    return;
}

/*!
    This function computes the modified solution for the well-balanced treatment of the
    gravitational source terms. The modified solution vector is given by
    \f{equation}{
      {\bf u}^* = \left[\begin{array}{c} \rho \varrho^{-1}\left(x,y\right) \\ \rho u \varrho^{-1}\left(x,y\right) \\ \rho v \varrho^{-1}\left(x,y\right) \\ e^* \end{array}\right]
    \f}
    where
    \f{equation}{
      e^* = \frac{p \varphi^{-1}\left(x,y\right)}{\gamma-1} + \frac{1}{2}\rho \varrho^{-1}\left(x,y\right) \left(u^2+v^2\right)
    \f}
    \f$\varrho\f$ and \f$\varphi\f$ are computed in #NavierStokes2DGravityField(). For flows without gravity, \f${\bf u}^* = {\bf u}\f$.

    References:
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, AIAA Journal, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580.
*/
extern "C" int gpuNavierStokes2DModifiedSolution(
    double  *uC,  /*!< Array to hold the computed modified solution */
    double  *u,   /*!< Solution vector array */
    int     d,    /*!< spatial dimension (not used) */
    void    *s,   /*!< Solver object of type #HyPar */
    void    *m,   /*!< MPI object of time #MPIVariables */
    double waqt   /*!< Current simulation time */
  )
{
    HyPar           *solver = (HyPar*)          s;
    NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;

    int     ghosts  = solver->ghosts;
    int     *dim    = solver->dim_local;
    int     ndims   = solver->ndims;

    int size = 1; for (int i=0; i <ndims; i++) size *= (dim[i] + 2*ghosts);
    int nblocks = (size - 1) / GPU_THREADS_PER_BLOCK + 1;
    double inv_gamma_m1 = 1.0 / (param->gamma-1.0);

    double cpu_time = 0.0;
    clock_t cpu_start, cpu_end;

    cpu_start = clock();
    NavierStokes2DModifiedSolution_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        size, param->gamma, inv_gamma_m1,
        param->gpu_grav_field_f, param->gpu_grav_field_g,
        u, uC);
    cudaDeviceSynchronize();
    cpu_end = clock();
    cpu_time += (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;

    return(0);
}

#endif
