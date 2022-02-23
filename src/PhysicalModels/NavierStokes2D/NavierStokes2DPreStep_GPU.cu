/*! @file NavierStokes2DPreStep_GPU.cu
    @brief Pre-step function for 2D Navier Stokes equations
    @author Youngdae Kim
*/

#include <time.h>
#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>
#include <matmult_native.h>
#include <basic_gpu.h>
#include <physicalmodels/navierstokes2d.h>
#include <hypar.h>

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel for gpuNavierStokes2DPreStep */
__global__
void NavierStokes2DPreStep_kernel(
    int size,
    double gamma,
    const double *u,
    double *fast_jac
)
{
    const int  JacSize = _MODEL_NVARS_*_MODEL_NVARS_;
    int        dir, p;
    double     *A;
    double     D[_MODEL_NVARS_*_MODEL_NVARS_],L[_MODEL_NVARS_*_MODEL_NVARS_],
               R[_MODEL_NVARS_*_MODEL_NVARS_],DL[_MODEL_NVARS_*_MODEL_NVARS_];

    p = threadIdx.x + blockDim.x * blockIdx.x;
    if (p < size) {
        dir = _XDIR_;
        A = (fast_jac + 2*JacSize*p + dir*JacSize);
        /* get the eigenvalues, and left & right eigenvectors */
        _NavierStokes2DEigenvalues_      ((u+_MODEL_NVARS_*p),D,gamma,dir);
        _NavierStokes2DLeftEigenvectors_ ((u+_MODEL_NVARS_*p),L,gamma,dir);
        _NavierStokes2DRightEigenvectors_((u+_MODEL_NVARS_*p),R,gamma,dir);
        /* remove the entropy modes (corresponding to eigenvalues u) */
        D[2*_MODEL_NVARS_+2] = D[3*_MODEL_NVARS_+3] = 0.0;
        /* assemble the Jacobian */
        MatMult4(_MODEL_NVARS_,DL,D,L );
        MatMult4(_MODEL_NVARS_,A,R,DL);

        dir = _YDIR_;
        A = (fast_jac + 2*JacSize*p + dir*JacSize);
        /* get the eigenvalues, and left & right eigenvectors */
        _NavierStokes2DEigenvalues_      ((u+_MODEL_NVARS_*p),D,gamma,dir)
        _NavierStokes2DLeftEigenvectors_ ((u+_MODEL_NVARS_*p),L,gamma,dir);
        _NavierStokes2DRightEigenvectors_((u+_MODEL_NVARS_*p),R,gamma,dir);
        /* remove the entropy modes (corresponding to eigenvalues v) */
        D[1*_MODEL_NVARS_+1] = D[3*_MODEL_NVARS_+3] = 0.0;
        /* assemble the Jacobian */
        MatMult4(_MODEL_NVARS_,DL,D,L );
        MatMult4(_MODEL_NVARS_,A,R,DL);
    }

    return;
}

/*! Pre-step function for the 2D Navier Stokes equations: This function
    is called at the beginning of each time step.
    + The solution at the beginning of the step is copied into #NavierStokes2D::solution
      for the linearized flux partitioning.
    + The linearized "fast" Jacobian (representing the acoustic modes) #NavierStokes2D::fast_jac
      is computed as:
      \f{equation}{
        A_f\left({\bf u}\right) = R\left({\bf u}\right)\Lambda_f\left({\bf u}\right)L\left({\bf u}\right)
      \f}
      where \f$R\f$ and \f$L\f$ are the matrices of right and left eigenvectors, and,
      \f{equation}{
        \Lambda_f = diag\left[0,0,u+c,u-c \right]
      \f}
      with \f$c\f$ being the speed of sound.
    \n\n

    Reference:
    + Ghosh, D., Constantinescu, E. M., Semi-Implicit Time Integration of Atmospheric Flows
      with Characteristic-Based Flux Partitioning, SIAM Journal on Scientific Computing,
      38 (3), 2016, A1848-A1875, http://dx.doi.org/10.1137/15M1044369.
*/
extern "C" int gpuNavierStokes2DPreStep(
    double  *u,   /*!< Solution vector */
    void    *s,   /*!< Solver object of type #HyPar */
    void    *m,   /*!< MPI object of type #MPIVariables */
    double  waqt  /*!< Current simulation time */
)
{
    HyPar             *solver    = (HyPar*)   s;
    NavierStokes2D    *param     = (NavierStokes2D*) solver->physics;
    int               *dim       = solver->dim_local;
    int               ghosts     = solver->ghosts;

    static int bounds[_MODEL_NDIMS_];
    _ArrayAddCopy1D_(dim,(2*ghosts),bounds,_MODEL_NDIMS_);
    int N_grid; _ArrayProduct1D_(bounds,_MODEL_NDIMS_,N_grid);
    int size = 2*N_grid*_MODEL_NVARS_*_MODEL_NVARS_;

    double cpu_time = 0.0;
    clock_t cpu_start, cpu_end;

    gpuMemcpy(solver->gpu_u, u, N_grid*_MODEL_NVARS_*sizeof(double), gpuMemcpyHostToDevice);
    gpuMemcpy(param->gpu_solution, u, N_grid*_MODEL_NVARS_*sizeof(double), gpuMemcpyHostToDevice);

    int nblocks = (N_grid - 1) / GPU_THREADS_PER_BLOCK + 1;
    cpu_start = clock();
    NavierStokes2DPreStep_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(N_grid, param->gamma, solver->gpu_u, param->gpu_fast_jac);
    cudaDeviceSynchronize();
    cpu_end = clock();
    cpu_time += (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;

    double *fast_jac_h = (double *)malloc(size*sizeof(double));
    gpuMemcpy(fast_jac_h, param->gpu_fast_jac, size*sizeof(double), gpuMemcpyDeviceToHost);

    double err = 0;
    for (int i = 0; i < size; i++) {
        if (fabs(fast_jac_h[i] - param->fast_jac[i]) > 1e-8)
            err += fabs(fast_jac_h[i] - param->fast_jac[i]);
    }
    free(fast_jac_h);

    return 0;
}

#endif
