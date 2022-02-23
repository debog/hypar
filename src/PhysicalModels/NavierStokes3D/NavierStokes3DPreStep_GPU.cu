/*! @file NavierStokes3DPreStep_GPU.cu
    @brief Pre-step function for 3D Navier Stokes equations
    @author Youngdae Kim
*/
#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>
#include <mathfunctions.h>
#include <matmult_native.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel for gpuNavierStokes3DPreStep() */
__global__
void gpuNavierStokes3DPreStep_kernel(
  int    npoints_grid,
  double gamma,
  const double * __restrict__ u,
  double       * __restrict__ fast_jac
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    double * __restrict__ A;
    double D[_MODEL_NVARS_*_MODEL_NVARS_],L[_MODEL_NVARS_*_MODEL_NVARS_],
           R[_MODEL_NVARS_*_MODEL_NVARS_],DL[_MODEL_NVARS_*_MODEL_NVARS_];
    int JacSize = _MODEL_NVARS_*_MODEL_NVARS_;
    int q = _MODEL_NVARS_*p;
    int dir = _XDIR_;

    A = (fast_jac + _MODEL_NDIMS_*JacSize*p + dir*JacSize);
    /* get the eigenvalues, and left & right eigenvectors */
    _NavierStokes3DEigenvalues_      ((u+q),_NavierStokes3D_stride_,D,gamma,dir);
    _NavierStokes3DLeftEigenvectors_ ((u+q),_NavierStokes3D_stride_,L,gamma,dir);
    _NavierStokes3DRightEigenvectors_((u+q),_NavierStokes3D_stride_,R,gamma,dir);
    /* remove the entropy modes (corresponding to eigenvalues u) */
    D[0] = D[12] = D[18] = 0.0;
    /* assemble the Jacobian */
    MatMult5(_MODEL_NVARS_,DL,D,L );
    MatMult5(_MODEL_NVARS_,A ,R,DL);

    dir = _YDIR_;
    A = (fast_jac + _MODEL_NDIMS_*JacSize*p + dir*JacSize);
    /* get the eigenvalues, and left & right eigenvectors */
    _NavierStokes3DEigenvalues_      ((u+q),_NavierStokes3D_stride_,D,gamma,dir)
    _NavierStokes3DLeftEigenvectors_ ((u+q),_NavierStokes3D_stride_,L,gamma,dir);
    _NavierStokes3DRightEigenvectors_((u+q),_NavierStokes3D_stride_,R,gamma,dir);
    /* remove the entropy modes (corresponding to eigenvalues v) */
    D[0] = D[6] = D[18] = 0.0;
    /* assemble the Jacobian */
    MatMult5(_MODEL_NVARS_,DL,D,L );
    MatMult5(_MODEL_NVARS_,A ,R,DL);

    dir = _ZDIR_;
    A = (fast_jac + _MODEL_NDIMS_*JacSize*p + dir*JacSize);
    /* get the eigenvalues, and left & right eigenvectors */
    _NavierStokes3DEigenvalues_      ((u+q),_NavierStokes3D_stride_,D,gamma,dir)
    _NavierStokes3DLeftEigenvectors_ ((u+q),_NavierStokes3D_stride_,L,gamma,dir);
    _NavierStokes3DRightEigenvectors_((u+q),_NavierStokes3D_stride_,R,gamma,dir);
    /* remove the entropy modes (corresponding to eigenvalues v) */
    D[0] = D[6] = D[12] = 0.0;
    /* assemble the Jacobian */
    MatMult5(_MODEL_NVARS_,DL,D,L );
    MatMult5(_MODEL_NVARS_,A ,R,DL);
  }
  return;
}

/*! Pre-step function for the 3D Navier Stokes equations: This function
    is called at the beginning of each time step.
    + The solution at the beginning of the step is copied into #NavierStokes3D::solution
      for the linearized flux partitioning.
    + The linearized "fast" Jacobian (representing the acoustic modes) #NavierStokes3D::fast_jac
      is computed as:
      \f{equation}{
        A_f\left({\bf u}\right) = R\left({\bf u}\right)\Lambda_f\left({\bf u}\right)L\left({\bf u}\right)
      \f}
      where \f$R\f$ and \f$L\f$ are the matrices of right and left eigenvectors, and,
      \f{equation}{
        \Lambda_f = diag\left[0,0,0,u+c,u-c \right]
      \f}
      with \f$c\f$ being the speed of sound.
    \n\n

    Reference:
    + Ghosh, D., Constantinescu, E. M., Semi-Implicit Time Integration of Atmospheric Flows
      with Characteristic-Based Flux Partitioning, SIAM Journal on Scientific Computing,
      38 (3), 2016, A1848-A1875, http://dx.doi.org/10.1137/15M1044369.
*/
extern "C" int gpuNavierStokes3DPreStep(
  double  *u,   /*!< Solution vector */
  void    *s,   /*!< Solver object of type #HyPar */
  void    *m,   /*!< MPI object of type #MPIVariables */
  double  waqt  /*!< Current simulation time */
)
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes3D    *param  = (NavierStokes3D*) solver->physics;

  int nblocks = (solver->npoints_local_wghosts-1)/GPU_THREADS_PER_BLOCK + 1;

  gpuArrayCopy1D(u,param->gpu_solution,(solver->npoints_local_wghosts)*_MODEL_NVARS_);
  gpuNavierStokes3DPreStep_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    solver->npoints_local_wghosts, param->gamma, u, param->gpu_fast_jac
  );
  cudaDeviceSynchronize();

  return 0;
}

#else

/*! Kernel for gpuNavierStokes3DPreStep() */
__global__
void gpuNavierStokes3DPreStep_kernel(
  int    npoints_grid,
  double gamma,
  const double * __restrict__ u,
  double       * __restrict__ fast_jac
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    double * __restrict__ A;
    double D[_MODEL_NVARS_*_MODEL_NVARS_],L[_MODEL_NVARS_*_MODEL_NVARS_],
           R[_MODEL_NVARS_*_MODEL_NVARS_],DL[_MODEL_NVARS_*_MODEL_NVARS_];
    int JacSize = _MODEL_NVARS_*_MODEL_NVARS_;
    int dir = _XDIR_;

    A = (fast_jac + _MODEL_NDIMS_*JacSize*p + dir*JacSize);
    /* get the eigenvalues, and left & right eigenvectors */
    _NavierStokes3DEigenvalues_      ((u+p),npoints_grid,D,gamma,dir);
    _NavierStokes3DLeftEigenvectors_ ((u+p),npoints_grid,L,gamma,dir);
    _NavierStokes3DRightEigenvectors_((u+p),npoints_grid,R,gamma,dir);
    /* remove the entropy modes (corresponding to eigenvalues u) */
    D[0] = D[12] = D[18] = 0.0;
    /* assemble the Jacobian */
    MatMult5(_MODEL_NVARS_,DL,D,L );
    MatMult5(_MODEL_NVARS_,A ,R,DL);

    dir = _YDIR_;
    A = (fast_jac + _MODEL_NDIMS_*JacSize*p + dir*JacSize);
    /* get the eigenvalues, and left & right eigenvectors */
    _NavierStokes3DEigenvalues_      ((u+p),npoints_grid,D,gamma,dir)
    _NavierStokes3DLeftEigenvectors_ ((u+p),npoints_grid,L,gamma,dir);
    _NavierStokes3DRightEigenvectors_((u+p),npoints_grid,R,gamma,dir);
    /* remove the entropy modes (corresponding to eigenvalues v) */
    D[0] = D[6] = D[18] = 0.0;
    /* assemble the Jacobian */
    MatMult5(_MODEL_NVARS_,DL,D,L );
    MatMult5(_MODEL_NVARS_,A ,R,DL);

    dir = _ZDIR_;
    A = (fast_jac + _MODEL_NDIMS_*JacSize*p + dir*JacSize);
    /* get the eigenvalues, and left & right eigenvectors */
    _NavierStokes3DEigenvalues_      ((u+p),npoints_grid,D,gamma,dir)
    _NavierStokes3DLeftEigenvectors_ ((u+p),npoints_grid,L,gamma,dir);
    _NavierStokes3DRightEigenvectors_((u+p),npoints_grid,R,gamma,dir);
    /* remove the entropy modes (corresponding to eigenvalues v) */
    D[0] = D[6] = D[12] = 0.0;
    /* assemble the Jacobian */
    MatMult5(_MODEL_NVARS_,DL,D,L );
    MatMult5(_MODEL_NVARS_,A ,R,DL);
  }
  return;
}

/*! Pre-step function for the 3D Navier Stokes equations: This function
    is called at the beginning of each time step.
    + The solution at the beginning of the step is copied into #NavierStokes3D::solution
      for the linearized flux partitioning.
    + The linearized "fast" Jacobian (representing the acoustic modes) #NavierStokes3D::fast_jac
      is computed as:
      \f{equation}{
        A_f\left({\bf u}\right) = R\left({\bf u}\right)\Lambda_f\left({\bf u}\right)L\left({\bf u}\right)
      \f}
      where \f$R\f$ and \f$L\f$ are the matrices of right and left eigenvectors, and,
      \f{equation}{
        \Lambda_f = diag\left[0,0,0,u+c,u-c \right]
      \f}
      with \f$c\f$ being the speed of sound.
    \n\n

    Reference:
    + Ghosh, D., Constantinescu, E. M., Semi-Implicit Time Integration of Atmospheric Flows
      with Characteristic-Based Flux Partitioning, SIAM Journal on Scientific Computing,
      38 (3), 2016, A1848-A1875, http://dx.doi.org/10.1137/15M1044369.
*/
extern "C" int gpuNavierStokes3DPreStep(
  double  *u,   /*!< Solution vector */
  void    *s,   /*!< Solver object of type #HyPar */
  void    *m,   /*!< MPI object of type #MPIVariables */
  double  waqt  /*!< Current simulation time */
)
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes3D    *param  = (NavierStokes3D*) solver->physics;

  int nblocks = (solver->npoints_local_wghosts-1)/GPU_THREADS_PER_BLOCK + 1;

  gpuArrayCopy1D(u,param->gpu_solution,(solver->npoints_local_wghosts)*_MODEL_NVARS_);
  gpuNavierStokes3DPreStep_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    solver->npoints_local_wghosts, param->gamma, u, param->gpu_fast_jac
  );
  cudaDeviceSynchronize();

  return 0;
}

#endif
