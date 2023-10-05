/*! @file NavierStokes2DUpwind_GPU.cu
    @author Youngdae Kim
    @brief Contains functions to compute the upwind flux at grid interfaces for the 2D Navier Stokes equations.
*/
#include <time.h>
#include <arrayfunctions_gpu.h>
#include <math_ops.h>
#include <basic_gpu.h>
#include <physicalmodels/navierstokes2d.h>
#include <hypar.h>

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel for gpuNavierStokes2DUpwindRusanov() */
__global__
void NavierStokes2DUpwindRusanov_kernel(
    int ngrid_points,
    int dir,
    int ghosts,
    double gamma,
    const int *dim,
    const double *grav_field_g,
    const double *fL,
    const double *fR,
    const double *uL,
    const double *uR,
    const double *u,
    double *fI
)
{
    int p = threadIdx.x + (blockDim.x * blockIdx.x);
    if (p < ngrid_points) {
        int bounds_inter[2], index_inter[2];
        bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[dir] += 1;
        _ArrayIndexnD_(_MODEL_NDIMS_,p,bounds_inter,index_inter,0);
        int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[dir]--;
        int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
        int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,ghosts,pL);
        int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,ghosts,pR);
        double udiff[_MODEL_NVARS_],uavg[_MODEL_NVARS_];

        /* Modified Rusanov's upwinding scheme */

        udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
        udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
        udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);
        udiff[3] = 0.5 * (uR[_MODEL_NVARS_*p+3] - uL[_MODEL_NVARS_*p+3]);

        _NavierStokes2DRoeAverage_(uavg,(u+_MODEL_NVARS_*pL),(u+_MODEL_NVARS_*pR),gamma);

        double c, vel[_MODEL_NDIMS_], rho,E,P;
        _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*pL),rho,vel[0],vel[1],E,P,gamma);
        c = sqrt(gamma*P/rho);
        double alphaL = c + absolute(vel[dir]);
        _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*pR),rho,vel[0],vel[1],E,P,gamma);
        c = sqrt(gamma*P/rho);
        double alphaR = c + absolute(vel[dir]);
        _NavierStokes2DGetFlowVar_(uavg,rho,vel[0],vel[1],E,P,gamma);
        c = sqrt(gamma*P/rho);
        double alphaavg = c + absolute(vel[dir]);

        double kappa  = max(grav_field_g[pL],grav_field_g[pR]);
        double alpha  = kappa*max3(alphaL,alphaR,alphaavg);

        fI[_MODEL_NVARS_*p+0] = 0.5*(fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0])-alpha*udiff[0];
        fI[_MODEL_NVARS_*p+1] = 0.5*(fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1])-alpha*udiff[1];
        fI[_MODEL_NVARS_*p+2] = 0.5*(fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2])-alpha*udiff[2];
        fI[_MODEL_NVARS_*p+3] = 0.5*(fL[_MODEL_NVARS_*p+3]+fR[_MODEL_NVARS_*p+3])-alpha*udiff[3];
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

    This upwinding scheme is modified for the balanced discretization of the 2D Navier Stokes equations when
    there is a non-zero gravitational force. See the reference below. For flows without any gravitational forces,
    it reduces to its original form.
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, AIAA Journal, http://dx.doi.org/10.2514/1.J054580.

*/
extern "C" int gpuNavierStokes2DUpwindRusanov(
    double *fI,
    double *fL,
    double *fR,
    double *uL,
    double *uR,
    double *u,
    int    dir,
    void   *s,
    double t
)
{
    HyPar          *solver = (HyPar *) s;
    NavierStokes2D *param  = (NavierStokes2D*) solver->physics;
    int            *dim    = solver->dim_local;

    int bounds_outer[2];
    bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[dir] = 1;
    int N_outer; _ArrayProduct1D_(bounds_outer,_MODEL_NDIMS_,N_outer);

    double cpu_time = 0.0;
    clock_t cpu_start, cpu_end;

    int ngrid_points = N_outer * (dim[dir]+1);
    int nblocks = (ngrid_points - 1) / GPU_THREADS_PER_BLOCK + 1;

    cpu_start = clock();
    NavierStokes2DUpwindRusanov_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        ngrid_points, dir, solver->ghosts, param->gamma, solver->dim_local,
        param->gpu_grav_field_g, fL, fR, uL, uR, u, fI);
    cudaDeviceSynchronize();
    cpu_end = clock();
    cpu_time += (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;

    return (0);
}

#endif
