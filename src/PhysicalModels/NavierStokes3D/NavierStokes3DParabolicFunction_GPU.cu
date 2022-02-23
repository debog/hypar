/*! @file NavierStokes3DParabolicFunction_GPU.cu
    @author Youngdae Kim
    @brief Compute the viscous terms for the 3D Navier Stokes equations
*/
#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <mpivars.h>
#include <hypar.h>

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel function */
__global__
void gpuNavierStokes3DParabolicFunction_Q_kernel(
    int npoints_grid,
    double gamma,
    const double * __restrict__ u,
    double * __restrict__ Q
)
{
    int p = blockDim.x * blockIdx.x + threadIdx.x;

    if (p < npoints_grid) {
        double energy, pressure;
        p *= _MODEL_NVARS_;
        _NavierStokes3DGetFlowVar_(
            (u+p),
            _NavierStokes3D_stride_,
            Q[p+0],
            Q[p+1],
            Q[p+2],
            Q[p+3],
            energy,
            pressure,
            gamma
        );
        Q[p+4] = gamma*pressure/Q[p+0]; /* temperature */
    }
    return;
}

/*! Kernel function */
__global__
void gpuNavierStokes3DParabolicFunction_QDeriv_kernel(
    int npoints_grid,
    int ghosts,
    const int * __restrict__ dim,
    const double * __restrict__ dxinv,
    double * __restrict__ QDerivX,
    double * __restrict__ QDerivY,
    double * __restrict__ QDerivZ
)
{
    int p = blockDim.x * blockIdx.x + threadIdx.x;

    if (p < npoints_grid) {
        int index[3];
        double dxinv_val, dyinv_val, dzinv_val;

        _ArrayIndexnD_(_MODEL_NDIMS_,p,dim,index,ghosts); p *= _MODEL_NVARS_;
        _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,dxinv,dxinv_val);
        _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,dxinv,dyinv_val);
        _GetCoordinate_(_ZDIR_,index[_ZDIR_],dim,ghosts,dxinv,dzinv_val);
        _ArrayScale1D_((QDerivX+p),dxinv_val,_MODEL_NVARS_);
        _ArrayScale1D_((QDerivY+p),dyinv_val,_MODEL_NVARS_);
        _ArrayScale1D_((QDerivZ+p),dzinv_val,_MODEL_NVARS_);
    }
    return;
}

/*! Kernel function */
__global__
void gpuNavierStokes3DParabolicFunction_Xdir_kernel(
    int npoints_grid,
    double two_third,
    double inv_gamma_m1,
    double inv_Re,
    double inv_Pr,
    const double * __restrict__ Q,
    const double * __restrict__ QDerivX,
    const double * __restrict__ QDerivY,
    const double * __restrict__ QDerivZ,
    double * __restrict__ FViscous
)
{
    int p = blockDim.x * blockIdx.x + threadIdx.x;

    if (p < npoints_grid) {
        double uvel, vvel, wvel, T, Tx,
               ux, uy, uz, vx, vy, wx, wz;

        p *= _MODEL_NVARS_;

        uvel = (Q+p)[1];
        vvel = (Q+p)[2];
        wvel = (Q+p)[3];
        T    = (Q+p)[4];
        Tx   = (QDerivX+p)[4];
        ux   = (QDerivX+p)[1];
        vx   = (QDerivX+p)[2];
        wx   = (QDerivX+p)[3];
        uy   = (QDerivY+p)[1];
        vy   = (QDerivY+p)[2];
        uz   = (QDerivZ+p)[1];
        wz   = (QDerivZ+p)[3];

        /* calculate viscosity coeff based on Sutherland's law */
        double mu = raiseto(T, 0.76);
        double tau_xx, tau_xy, tau_xz, qx;
        tau_xx = two_third * (mu*inv_Re) * (2*ux - vy - wz);
        tau_xy = (mu*inv_Re) * (uy + vx);
        tau_xz = (mu*inv_Re) * (uz + wx);
        qx     = ( mu*inv_Re * inv_gamma_m1 * inv_Pr ) * Tx;

        (FViscous+p)[0] = 0.0;
        (FViscous+p)[1] = tau_xx;
        (FViscous+p)[2] = tau_xy;
        (FViscous+p)[3] = tau_xz;
        (FViscous+p)[4] = uvel*tau_xx + vvel*tau_xy + wvel*tau_xz + qx;
    }
    return;
}

/*! Kernel function */
__global__
void gpuNavierStokes3DParabolicFunction_Ydir_kernel(
    int npoints_grid,
    double two_third,
    double inv_gamma_m1,
    double inv_Re,
    double inv_Pr,
    const double * __restrict__ Q,
    const double * __restrict__ QDerivX,
    const double * __restrict__ QDerivY,
    const double * __restrict__ QDerivZ,
    double * __restrict__ FViscous
)
{
    int p = blockDim.x * blockIdx.x + threadIdx.x;
    if (p < npoints_grid) {
        double uvel, vvel, wvel, T, Ty,
               ux, uy, vx, vy, vz, wy, wz;

        p *= _MODEL_NVARS_;

        uvel = (Q+p)[1];
        vvel = (Q+p)[2];
        wvel = (Q+p)[3];
        T    = (Q+p)[4];
        Ty   = (QDerivY+p)[4];
        ux   = (QDerivX+p)[1];
        vx   = (QDerivX+p)[2];
        uy   = (QDerivY+p)[1];
        vy   = (QDerivY+p)[2];
        wy   = (QDerivY+p)[3];
        vz   = (QDerivZ+p)[2];
        wz   = (QDerivZ+p)[3];

        /* calculate viscosity coeff based on Sutherland's law */
        double mu = raiseto(T, 0.76);

        double tau_yx, tau_yy, tau_yz, qy;
        tau_yx = (mu*inv_Re) * (uy + vx);
        tau_yy = two_third * (mu*inv_Re) * (-ux + 2*vy - wz);
        tau_yz = (mu*inv_Re) * (vz + wy);
        qy     = ( mu*inv_Re * inv_gamma_m1 * inv_Pr ) * Ty;

        (FViscous+p)[0] = 0.0;
        (FViscous+p)[1] = tau_yx;
        (FViscous+p)[2] = tau_yy;
        (FViscous+p)[3] = tau_yz;
        (FViscous+p)[4] = uvel*tau_yx + vvel*tau_yy + wvel*tau_yz + qy;
    }
    return;
}

/*! Kernel function */
__global__
void gpuNavierStokes3DParabolicFunction_Zdir_kernel(
    int npoints_grid,
    double two_third,
    double inv_gamma_m1,
    double inv_Re,
    double inv_Pr,
    const double * __restrict__ Q,
    const double * __restrict__ QDerivX,
    const double * __restrict__ QDerivY,
    const double * __restrict__ QDerivZ,
    double * __restrict__ FViscous
)
{
    int p = blockDim.x * blockIdx.x + threadIdx.x;
    if (p < npoints_grid) {
        double uvel, vvel, wvel, T, Tz,
               ux, uz, vy, vz, wx, wy, wz;

        p *= _MODEL_NVARS_;

        uvel = (Q+p)[1];
        vvel = (Q+p)[2];
        wvel = (Q+p)[3];
        T    = (Q+p)[4];
        Tz   = (QDerivZ+p)[4];
        ux   = (QDerivX+p)[1];
        wx   = (QDerivX+p)[3];
        vy   = (QDerivY+p)[2];
        wy   = (QDerivY+p)[3];
        uz   = (QDerivZ+p)[1];
        vz   = (QDerivZ+p)[2];
        wz   = (QDerivZ+p)[3];

        /* calculate viscosity coeff based on Sutherland's law */
        double mu = raiseto(T,0.76);

        double tau_zx, tau_zy, tau_zz, qz;
        tau_zx = (mu*inv_Re) * (uz + wx);
        tau_zy = (mu*inv_Re) * (vz + wy);
        tau_zz = two_third * (mu*inv_Re) * (-ux - vy + 2*wz);
        qz     = ( mu*inv_Re * inv_gamma_m1 * inv_Pr ) * Tz;

        (FViscous+p)[0] = 0.0;
        (FViscous+p)[1] = tau_zx;
        (FViscous+p)[2] = tau_zy;
        (FViscous+p)[3] = tau_zz;
        (FViscous+p)[4] = uvel*tau_zx + vvel*tau_zy + wvel*tau_zz + qz;
    }
    return;
}

/*! Kernel function */
__global__
void gpuNavierStokes3DParabolicFunction_par_kernel(
    int npoints_grid,
    int ghosts,
    int dir,
    const int * __restrict__ dim,
    const double * __restrict__ dxinv,
    const double * __restrict__ FDeriv,
    double * __restrict__ par
)
{
    int p = blockDim.x * blockIdx.x + threadIdx.x;

    if (p < npoints_grid) {
        int index[3];
        double dxinv_val;

        _ArrayIndexnD_(_MODEL_NDIMS_,p,dim,index,0);
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p); p *= _MODEL_NVARS_;
        _GetCoordinate_(dir,index[dir],dim,ghosts,dxinv,dxinv_val);

        #pragma unroll
        for (int v=0; v<_MODEL_NVARS_; v++) (par+p)[v] += (dxinv_val * (FDeriv+p)[v] );
    }
    return;
}

/*!
    Compute the viscous terms in the 3D Navier Stokes equations: this function computes
    the following:
    \f{equation}{
      \frac {\partial} {\partial x} \left[\begin{array}{c} 0 \\ \tau_{xx} \\ \tau_{yx} \\ \tau_{zx} \\ u \tau_{xx} + v \tau_{yx} + w \tau_{zx} - q_x \end{array}\right]
      + \frac {\partial} {\partial y} \left[\begin{array}{c} 0 \\ \tau_{xy} \\ \tau_{yy} \\ \tau_{zy} \\ u \tau_{xy} + v \tau_{yy} + w \tau_{zy} - q_y \end{array}\right]
      + \frac {\partial} {\partial z} \left[\begin{array}{c} 0 \\ \tau_{xz} \\ \tau_{yz} \\ \tau_{zz} \\ u \tau_{xz} + v \tau_{yz} + w \tau_{zz} - q_z \end{array}\right]
    \f}
    where
    \f{align}{
      \tau_{xx} &= \frac{2}{3}\left(\frac{\mu}{Re}\right)\left(2\frac{\partial u}{\partial x} - \frac{\partial v}{\partial y} - \frac{\partial w}{\partial z}\right),\\
      \tau_{xy} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right),\\
      \tau_{xz} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial z} + \frac{\partial w}{\partial x}\right),\\
      \tau_{yx} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right),\\
      \tau_{yy} &= \frac{2}{3}\left(\frac{\mu}{Re}\right)\left(-\frac{\partial u}{\partial x} +2\frac{\partial v}{\partial y} - \frac{\partial w}{\partial z}\right),\\
      \tau_{yz} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial w}{\partial y} + \frac{\partial v}{\partial z}\right),\\
      \tau_{zx} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial z} + \frac{\partial w}{\partial x}\right),\\
      \tau_{zy} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial v}{\partial z} + \frac{\partial w}{\partial y}\right),\\
      \tau_{zz} &= \frac{2}{3}\left(\frac{\mu}{Re}\right)\left(-\frac{\partial u}{\partial x} - \frac{\partial v}{\partial y} + 2\frac{\partial v}{\partial y}\right),\\
      q_x &= -\frac{mu}{\left(\gamma-1\right)Re Pr}\frac{\partial T}{\partial x}, \\
      q_y &= -\frac{mu}{\left(\gamma-1\right)Re Pr}\frac{\partial T}{\partial y}, \\
      q_z &= -\frac{mu}{\left(\gamma-1\right)Re Pr}\frac{\partial T}{\partial z}
    \f}
    and the temperature is \f$T = \gamma p/\rho\f$. \f$Re\f$ and \f$Pr\f$ are the Reynolds and Prandtl numbers, respectively. Note that this function
    computes the entire parabolic term, and thus bypasses HyPar's parabolic function calculation interfaces. NavierStokes3DInitialize() assigns this
    function to #HyPar::ParabolicFunction.
    \n\n
    Reference:
    + Tannehill, Anderson and Pletcher, Computational Fluid Mechanics and Heat Transfer,
      Chapter 5, Section 5.1.7 (However, the non-dimensional velocity and the Reynolds
      number is based on speed of sound, instead of the freestream velocity).
*/
extern "C" int gpuNavierStokes3DParabolicFunction(
    double * __restrict__ par,  /*!< Array to hold the computed viscous terms */
    double * __restrict__ u,    /*!< Solution vector array */
    void   * __restrict__ s,    /*!< Solver object of type #HyPar */
    void   * __restrict__ m,    /*!< MPI object of type #MPIVariables */
    double t                    /*!< Current simulation time */
)
{
  HyPar           *solver   = (HyPar*) s;
  MPIVariables    *mpi      = (MPIVariables*) m;
  NavierStokes3D  *physics  = (NavierStokes3D*) solver->physics;

  int ghosts = solver->ghosts;
  int *dim   = solver->gpu_dim_local;
  int size   = solver->npoints_local_wghosts;

  gpuMemset(par, 0, size*_MODEL_NVARS_*sizeof(double));
  if (physics->Re <= 0) return (0); /* inviscid flow */
  solver->count_par++;

  static double two_third    = 2.0/3.0;
  double        inv_gamma_m1 = 1.0 / (physics->gamma-1.0);
  double        inv_Re       = 1.0 / physics->Re;
  double        inv_Pr       = 1.0 / physics->Pr;

  double *Q       = physics->gpu_Q;
  double *QDerivX = physics->gpu_QDerivX;
  double *QDerivY = physics->gpu_QDerivY;
  double *QDerivZ = physics->gpu_QDerivZ;
  double *dxinv   = solver->gpu_dxinv;

  int nblocks = (size-1)/GPU_THREADS_PER_BLOCK + 1;
  int nblocks_par = (solver->npoints_local-1)/GPU_THREADS_PER_BLOCK + 1;

#if defined(GPU_STAT)
  cudaEvent_t startEvent, stopEvent;
  float milliseconds;

  checkCuda( cudaEventCreate(&startEvent) );
  checkCuda( cudaEventCreate(&stopEvent) );
  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  gpuNavierStokes3DParabolicFunction_Q_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      size, physics->gamma, u, Q
  );
#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
#endif

  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  int memory_accessed = size*_MODEL_NVARS_*sizeof(double);
  printf("%-50s GPU time (secs) = %.6f bandwidth (GB/s) = %6.2f\n",
          "NavierStokes3DParabolicFunction_Q", milliseconds*1e-3,
          (1e-6*(memory_accessed)/milliseconds));
#endif

  solver->FirstDerivativePar(QDerivX,Q,_XDIR_,1,solver,mpi);
  solver->FirstDerivativePar(QDerivY,Q,_YDIR_,1,solver,mpi);
  solver->FirstDerivativePar(QDerivZ,Q,_ZDIR_,1,solver,mpi);

  gpuMPIExchangeBoundariesnD(_MODEL_NDIMS_,_MODEL_NVARS_,solver->gpu_dim_local,
                             solver->ghosts,mpi,QDerivX);
  gpuMPIExchangeBoundariesnD(_MODEL_NDIMS_,_MODEL_NVARS_,solver->gpu_dim_local,
                             solver->ghosts,mpi,QDerivY);
  gpuMPIExchangeBoundariesnD(_MODEL_NDIMS_,_MODEL_NVARS_,solver->gpu_dim_local,
                             solver->ghosts,mpi,QDerivY);

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  gpuNavierStokes3DParabolicFunction_QDeriv_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      size, ghosts, dim, dxinv, QDerivX, QDerivY, QDerivZ
  );
#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
#endif

  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  memory_accessed = (3*size*_MODEL_NVARS_ + solver->size_x)*sizeof(double);
  printf("%-50s GPU time (secs) = %.6f bandwidth (GB/s) = %6.2f\n",
          "NavierStokes3DParabolicFunction_QDeriv", milliseconds*1e-3,
          (1e-6*(memory_accessed)/milliseconds));
#endif

  double *FViscous = physics->gpu_FViscous;
  double *FDeriv   = physics->gpu_FDeriv;
  gpuMemset(FViscous, 0, size*_MODEL_NVARS_*sizeof(double));
  gpuMemset(FDeriv, 0, size*_MODEL_NVARS_*sizeof(double));

  /* Along X */
#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  gpuNavierStokes3DParabolicFunction_Xdir_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      size, two_third, inv_gamma_m1, inv_Re, inv_Pr, Q,
      QDerivX, QDerivY, QDerivZ, FViscous
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  memory_accessed = 17*size*sizeof(double);
  printf("%-50s GPU time (secs) = %.6f bandwidth (GB/s) = %6.2f\n",
          "NavierStokes3DParabolicFunction_Xdir", milliseconds*1e-3,
          (1e-6*(memory_accessed)/milliseconds));
#endif

  solver->FirstDerivativePar(FDeriv,FViscous,_XDIR_,1,solver,mpi);

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  gpuNavierStokes3DParabolicFunction_par_kernel<<<nblocks_par, GPU_THREADS_PER_BLOCK>>>(
      solver->npoints_local, ghosts, _XDIR_, dim, dxinv, FDeriv, par
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  memory_accessed = (solver->npoints_local + 2*solver->npoints_local*_MODEL_NVARS_)*sizeof(double);
  printf("%-50s GPU time (secs) = %.6f bandwidth (GB/s) = %6.2f\n",
          "NavierStokes3DParabolicFunction_par", milliseconds*1e-3,
          (1e-6*(memory_accessed)/milliseconds));
#endif

  /* Along Y */
#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  gpuNavierStokes3DParabolicFunction_Ydir_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      size, two_third, inv_gamma_m1, inv_Re, inv_Pr, Q,
      QDerivX, QDerivY, QDerivZ, FViscous
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );
  //cudaDeviceSynchronize();
  memory_accessed = 17*size*sizeof(double);
  printf("%-50s GPU time (secs) = %.6f bandwidth (GB/s) = %6.2f\n",
          "NavierStokes3DParabolicFunction_Ydir", milliseconds*1e-3,
          (1e-6*(memory_accessed)/milliseconds));
#endif

  solver->FirstDerivativePar(FDeriv,FViscous,_YDIR_,1,solver,mpi);
  gpuNavierStokes3DParabolicFunction_par_kernel<<<nblocks_par, GPU_THREADS_PER_BLOCK>>>(
      solver->npoints_local, ghosts, _YDIR_, dim, dxinv, FDeriv, par
  );
  cudaDeviceSynchronize();

  /* Along Z */
#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  gpuNavierStokes3DParabolicFunction_Zdir_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      size, two_third, inv_gamma_m1, inv_Re, inv_Pr, Q,
      QDerivX, QDerivY, QDerivZ, FViscous
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  printf("%-50s GPU time (secs) = %.6f bandwidth (GB/s) = %6.2f\n",
          "NavierStokes3DParabolicFunction_Zdir", milliseconds*1e-3,
          (1e-6*(memory_accessed)/milliseconds));
#endif

  solver->FirstDerivativePar(FDeriv,FViscous,_ZDIR_,1,solver,mpi);
  gpuNavierStokes3DParabolicFunction_par_kernel<<<nblocks_par, GPU_THREADS_PER_BLOCK>>>(
      solver->npoints_local, ghosts, _ZDIR_, dim, dxinv, FDeriv, par
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventDestroy(startEvent) );
  checkCuda( cudaEventDestroy(stopEvent) );
#endif

  if (solver->flag_ib) gpuArrayBlockMultiply(par,solver->gpu_iblank,size,_MODEL_NVARS_);

  return (0);
}

#else

/*! Kernel function */
__global__
void gpuNavierStokes3DParabolicFunction_Q_kernel(
    int npoints_grid,
    double gamma,
    const double * __restrict__ u,
    double * __restrict__ Q
)
{
    int p = blockDim.x * blockIdx.x + threadIdx.x;

    if (p < npoints_grid) {
        double energy, pressure;
        _NavierStokes3DGetFlowVar_(
            (u+p),
            npoints_grid,
            Q[p],
            Q[p+  npoints_grid],
            Q[p+2*npoints_grid],
            Q[p+3*npoints_grid],
            energy,
            pressure,
            gamma
        );
        Q[p+4*npoints_grid] = gamma*pressure/Q[p+0]; /* temperature */
    }
    return;
}

/*! Kernel function */
__global__
void gpuNavierStokes3DParabolicFunction_QDeriv_kernel(
    int npoints_grid,
    int ghosts,
    const int * __restrict__ dim,
    const double * __restrict__ dxinv,
    double * __restrict__ QDerivX,
    double * __restrict__ QDerivY,
    double * __restrict__ QDerivZ
)
{
    int p = blockDim.x * blockIdx.x + threadIdx.x;

    if (p < npoints_grid) {
        int index[_MODEL_NDIMS_];
        double dxinv_val, dyinv_val, dzinv_val;

        _ArrayIndexnD_(_MODEL_NDIMS_,p,dim,index,ghosts);
        _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,dxinv,dxinv_val);
        _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,dxinv,dyinv_val);
        _GetCoordinate_(_ZDIR_,index[_ZDIR_],dim,ghosts,dxinv,dzinv_val);

        #pragma unroll
        for (int j=0; j<_MODEL_NVARS_; j++) {
          QDerivX[p] *= dxinv_val;
          QDerivY[p] *= dyinv_val;
          QDerivZ[p] *= dzinv_val;
          p += npoints_grid;
        }
    }
    return;
}

/*! Kernel function */
__global__
void gpuNavierStokes3DParabolicFunction_Xdir_kernel(
    int npoints_grid,
    double two_third,
    double inv_gamma_m1,
    double inv_Re,
    double inv_Pr,
    const double * __restrict__ Q,
    const double * __restrict__ QDerivX,
    const double * __restrict__ QDerivY,
    const double * __restrict__ QDerivZ,
    double * __restrict__ FViscous
)
{
    int p = blockDim.x * blockIdx.x + threadIdx.x;

    if (p < npoints_grid) {
        double uvel, vvel, wvel, T, Tx,
               ux, uy, uz, vx, vy, wx, wz;

        uvel = (Q+p+  npoints_grid)[0];
        vvel = (Q+p+2*npoints_grid)[0];
        wvel = (Q+p+3*npoints_grid)[0];
        T    = (Q+p+4*npoints_grid)[0];
        Tx   = (QDerivX+p+4*npoints_grid)[0];
        ux   = (QDerivX+p+  npoints_grid)[0];
        vx   = (QDerivX+p+2*npoints_grid)[0];
        wx   = (QDerivX+p+3*npoints_grid)[0];
        uy   = (QDerivY+p+  npoints_grid)[0];
        vy   = (QDerivY+p+2*npoints_grid)[0];
        uz   = (QDerivZ+p+  npoints_grid)[0];
        wz   = (QDerivZ+p+3*npoints_grid)[0];

        /* calculate viscosity coeff based on Sutherland's law */
        double mu = raiseto(T, 0.76);
        double tau_xx, tau_xy, tau_xz, qx;
        tau_xx = two_third * (mu*inv_Re) * (2*ux - vy - wz);
        tau_xy = (mu*inv_Re) * (uy + vx);
        tau_xz = (mu*inv_Re) * (uz + wx);
        qx     = ( mu*inv_Re * inv_gamma_m1 * inv_Pr ) * Tx;

        (FViscous+p)[0] = 0.0;
        (FViscous+p+  npoints_grid)[0] = tau_xx;
        (FViscous+p+2*npoints_grid)[0] = tau_xy;
        (FViscous+p+3*npoints_grid)[0] = tau_xz;
        (FViscous+p+4*npoints_grid)[0] = uvel*tau_xx + vvel*tau_xy + wvel*tau_xz + qx;
    }
    return;
}

/*! Kernel function */
__global__
void gpuNavierStokes3DParabolicFunction_Ydir_kernel(
    int npoints_grid,
    double two_third,
    double inv_gamma_m1,
    double inv_Re,
    double inv_Pr,
    const double * __restrict__ Q,
    const double * __restrict__ QDerivX,
    const double * __restrict__ QDerivY,
    const double * __restrict__ QDerivZ,
    double * __restrict__ FViscous
)
{
    int p = blockDim.x * blockIdx.x + threadIdx.x;
    if (p < npoints_grid) {
        double uvel, vvel, wvel, T, Ty,
               ux, uy, vx, vy, vz, wy, wz;

        uvel = (Q+p+  npoints_grid)[0];
        vvel = (Q+p+2*npoints_grid)[0];
        wvel = (Q+p+3*npoints_grid)[0];
        T    = (Q+p+4*npoints_grid)[0];
        Ty   = (QDerivY+p+4*npoints_grid)[0];
        ux   = (QDerivX+p+  npoints_grid)[0];
        vx   = (QDerivX+p+2*npoints_grid)[0];
        uy   = (QDerivY+p+  npoints_grid)[0];
        vy   = (QDerivY+p+2*npoints_grid)[0];
        wy   = (QDerivY+p+3*npoints_grid)[0];
        vz   = (QDerivZ+p+2*npoints_grid)[0];
        wz   = (QDerivZ+p+3*npoints_grid)[0];

        /* calculate viscosity coeff based on Sutherland's law */
        double mu = raiseto(T, 0.76);

        double tau_yx, tau_yy, tau_yz, qy;
        tau_yx = (mu*inv_Re) * (uy + vx);
        tau_yy = two_third * (mu*inv_Re) * (-ux + 2*vy - wz);
        tau_yz = (mu*inv_Re) * (vz + wy);
        qy     = ( mu*inv_Re * inv_gamma_m1 * inv_Pr ) * Ty;

        (FViscous+p)[0] = 0.0;
        (FViscous+p+  npoints_grid)[0] = tau_yx;
        (FViscous+p+2*npoints_grid)[0] = tau_yy;
        (FViscous+p+3*npoints_grid)[0] = tau_yz;
        (FViscous+p+4*npoints_grid)[0] = uvel*tau_yx + vvel*tau_yy + wvel*tau_yz + qy;
    }
    return;
}

/*! Kernel function */
__global__
void gpuNavierStokes3DParabolicFunction_Zdir_kernel(
    int npoints_grid,
    double two_third,
    double inv_gamma_m1,
    double inv_Re,
    double inv_Pr,
    const double * __restrict__ Q,
    const double * __restrict__ QDerivX,
    const double * __restrict__ QDerivY,
    const double * __restrict__ QDerivZ,
    double * __restrict__ FViscous
)
{
    int p = blockDim.x * blockIdx.x + threadIdx.x;
    if (p < npoints_grid) {
        double uvel, vvel, wvel, T, Tz,
               ux, uz, vy, vz, wx, wy, wz;

        uvel = (Q+p+  npoints_grid)[0];
        vvel = (Q+p+2*npoints_grid)[0];
        wvel = (Q+p+3*npoints_grid)[0];
        T    = (Q+p+4*npoints_grid)[0];
        Tz   = (QDerivZ+p+4*npoints_grid)[0];
        ux   = (QDerivX+p+  npoints_grid)[0];
        wx   = (QDerivX+p+3*npoints_grid)[0];
        vy   = (QDerivY+p+2*npoints_grid)[0];
        wy   = (QDerivY+p+3*npoints_grid)[0];
        uz   = (QDerivZ+p+  npoints_grid)[0];
        vz   = (QDerivZ+p+2*npoints_grid)[0];
        wz   = (QDerivZ+p+3*npoints_grid)[0];

        /* calculate viscosity coeff based on Sutherland's law */
        double mu = raiseto(T,0.76);

        double tau_zx, tau_zy, tau_zz, qz;
        tau_zx = (mu*inv_Re) * (uz + wx);
        tau_zy = (mu*inv_Re) * (vz + wy);
        tau_zz = two_third * (mu*inv_Re) * (-ux - vy + 2*wz);
        qz     = ( mu*inv_Re * inv_gamma_m1 * inv_Pr ) * Tz;

        (FViscous+p)[0] = 0.0;
        (FViscous+p+  npoints_grid)[0] = tau_zx;
        (FViscous+p+2*npoints_grid)[0] = tau_zy;
        (FViscous+p+3*npoints_grid)[0] = tau_zz;
        (FViscous+p+4*npoints_grid)[0] = uvel*tau_zx + vvel*tau_zy + wvel*tau_zz + qz;
    }
    return;
}

/*! Kernel function */
__global__
void gpuNavierStokes3DParabolicFunction_par_kernel(
    int npoints_grid,
    int npoints_local_wghosts,
    int ghosts,
    int dir,
    const int * __restrict__ dim,
    const double * __restrict__ dxinv,
    const double * __restrict__ FDeriv,
    double * __restrict__ par
)
{
    int p = blockDim.x * blockIdx.x + threadIdx.x;

    if (p < npoints_grid) {
        int index[3];
        double dxinv_val;

        _ArrayIndexnD_(_MODEL_NDIMS_,p,dim,index,0);
        _ArrayIndex1D_(_MODEL_NDIMS_,dim,index,ghosts,p);
        _GetCoordinate_(dir,index[dir],dim,ghosts,dxinv,dxinv_val);

        #pragma unroll
        for (int v=0; v<_MODEL_NVARS_; v++) {
            (par+p+v*npoints_local_wghosts)[0] += (dxinv_val * (FDeriv+p+v*npoints_local_wghosts)[0]);
        }
    }
    return;
}

/*!
    Compute the viscous terms in the 3D Navier Stokes equations: this function computes
    the following:
    \f{equation}{
      \frac {\partial} {\partial x} \left[\begin{array}{c} 0 \\ \tau_{xx} \\ \tau_{yx} \\ \tau_{zx} \\ u \tau_{xx} + v \tau_{yx} + w \tau_{zx} - q_x \end{array}\right]
      + \frac {\partial} {\partial y} \left[\begin{array}{c} 0 \\ \tau_{xy} \\ \tau_{yy} \\ \tau_{zy} \\ u \tau_{xy} + v \tau_{yy} + w \tau_{zy} - q_y \end{array}\right]
      + \frac {\partial} {\partial z} \left[\begin{array}{c} 0 \\ \tau_{xz} \\ \tau_{yz} \\ \tau_{zz} \\ u \tau_{xz} + v \tau_{yz} + w \tau_{zz} - q_z \end{array}\right]
    \f}
    where
    \f{align}{
      \tau_{xx} &= \frac{2}{3}\left(\frac{\mu}{Re}\right)\left(2\frac{\partial u}{\partial x} - \frac{\partial v}{\partial y} - \frac{\partial w}{\partial z}\right),\\
      \tau_{xy} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right),\\
      \tau_{xz} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial z} + \frac{\partial w}{\partial x}\right),\\
      \tau_{yx} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right),\\
      \tau_{yy} &= \frac{2}{3}\left(\frac{\mu}{Re}\right)\left(-\frac{\partial u}{\partial x} +2\frac{\partial v}{\partial y} - \frac{\partial w}{\partial z}\right),\\
      \tau_{yz} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial w}{\partial y} + \frac{\partial v}{\partial z}\right),\\
      \tau_{zx} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial z} + \frac{\partial w}{\partial x}\right),\\
      \tau_{zy} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial v}{\partial z} + \frac{\partial w}{\partial y}\right),\\
      \tau_{zz} &= \frac{2}{3}\left(\frac{\mu}{Re}\right)\left(-\frac{\partial u}{\partial x} - \frac{\partial v}{\partial y} + 2\frac{\partial v}{\partial y}\right),\\
      q_x &= -\frac{mu}{\left(\gamma-1\right)Re Pr}\frac{\partial T}{\partial x}, \\
      q_y &= -\frac{mu}{\left(\gamma-1\right)Re Pr}\frac{\partial T}{\partial y}, \\
      q_z &= -\frac{mu}{\left(\gamma-1\right)Re Pr}\frac{\partial T}{\partial z}
    \f}
    and the temperature is \f$T = \gamma p/\rho\f$. \f$Re\f$ and \f$Pr\f$ are the Reynolds and Prandtl numbers, respectively. Note that this function
    computes the entire parabolic term, and thus bypasses HyPar's parabolic function calculation interfaces. NavierStokes3DInitialize() assigns this
    function to #HyPar::ParabolicFunction.
    \n\n
    Reference:
    + Tannehill, Anderson and Pletcher, Computational Fluid Mechanics and Heat Transfer,
      Chapter 5, Section 5.1.7 (However, the non-dimensional velocity and the Reynolds
      number is based on speed of sound, instead of the freestream velocity).
*/
extern "C" int gpuNavierStokes3DParabolicFunction(
    double * __restrict__ par,  /*!< Array to hold the computed viscous terms */
    double * __restrict__ u,    /*!< Solution vector array */
    void   * __restrict__ s,    /*!< Solver object of type #HyPar */
    void   * __restrict__ m,    /*!< MPI object of type #MPIVariables */
    double t                    /*!< Current simulation time */
)
{
  HyPar           *solver   = (HyPar*) s;
  MPIVariables    *mpi      = (MPIVariables*) m;
  NavierStokes3D  *physics  = (NavierStokes3D*) solver->physics;

  int ghosts = solver->ghosts;
  int *dim   = solver->gpu_dim_local;
  int size   = solver->npoints_local_wghosts;

  gpuMemset(par, 0, size*_MODEL_NVARS_*sizeof(double));
  if (physics->Re <= 0) return (0); /* inviscid flow */
  solver->count_par++;

  static double two_third    = 2.0/3.0;
  double        inv_gamma_m1 = 1.0 / (physics->gamma-1.0);
  double        inv_Re       = 1.0 / physics->Re;
  double        inv_Pr       = 1.0 / physics->Pr;

  double *Q       = physics->gpu_Q;
  double *QDerivX = physics->gpu_QDerivX;
  double *QDerivY = physics->gpu_QDerivY;
  double *QDerivZ = physics->gpu_QDerivZ;
  double *dxinv   = solver->gpu_dxinv;

  int nblocks = (size-1)/GPU_THREADS_PER_BLOCK + 1;
  int nblocks_par = (solver->npoints_local-1)/GPU_THREADS_PER_BLOCK + 1;

#if defined(GPU_STAT)
  cudaEvent_t startEvent, stopEvent;
  float milliseconds;

  checkCuda( cudaEventCreate(&startEvent) );
  checkCuda( cudaEventCreate(&stopEvent) );

  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  gpuNavierStokes3DParabolicFunction_Q_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      size, physics->gamma, u, Q
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  int memory_accessed = size*_MODEL_NVARS_*sizeof(double);
  printf("%-50s GPU time (secs) = %.6f bandwidth (GB/s) = %6.2f\n",
          "NavierStokes3DParabolicFunction_Q", milliseconds*1e-3,
          (1e-6*(memory_accessed)/milliseconds));
#endif

  solver->FirstDerivativePar(QDerivX,Q,_XDIR_,1,solver,mpi);
  solver->FirstDerivativePar(QDerivY,Q,_YDIR_,1,solver,mpi);
  solver->FirstDerivativePar(QDerivZ,Q,_ZDIR_,1,solver,mpi);

  gpuMPIExchangeBoundariesnD(_MODEL_NDIMS_,_MODEL_NVARS_,solver->gpu_dim_local,
                             solver->ghosts,mpi,QDerivX);
  gpuMPIExchangeBoundariesnD(_MODEL_NDIMS_,_MODEL_NVARS_,solver->gpu_dim_local,
                             solver->ghosts,mpi,QDerivY);
  gpuMPIExchangeBoundariesnD(_MODEL_NDIMS_,_MODEL_NVARS_,solver->gpu_dim_local,
                              solver->ghosts,mpi,QDerivY);

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  gpuNavierStokes3DParabolicFunction_QDeriv_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      size, ghosts, dim, dxinv, QDerivX, QDerivY, QDerivZ
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  memory_accessed = (3*size*_MODEL_NVARS_ + solver->size_x)*sizeof(double);
  printf("%-50s GPU time (secs) = %.6f bandwidth (GB/s) = %6.2f\n",
          "NavierStokes3DParabolicFunction_QDeriv", milliseconds*1e-3,
          (1e-6*(memory_accessed)/milliseconds));
#endif

  double *FViscous = physics->gpu_FViscous;
  double *FDeriv   = physics->gpu_FDeriv;
  gpuMemset(FViscous, 0, size*_MODEL_NVARS_*sizeof(double));
  gpuMemset(FDeriv, 0, size*_MODEL_NVARS_*sizeof(double));

  /* Along X */

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  gpuNavierStokes3DParabolicFunction_Xdir_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      size, two_third, inv_gamma_m1, inv_Re, inv_Pr, Q,
      QDerivX, QDerivY, QDerivZ, FViscous
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  memory_accessed = 17*size*sizeof(double);
  printf("%-50s GPU time (secs) = %.6f bandwidth (GB/s) = %6.2f\n",
          "NavierStokes3DParabolicFunction_Xdir", milliseconds*1e-3,
          (1e-6*(memory_accessed)/milliseconds));
#endif

  solver->FirstDerivativePar(FDeriv,FViscous,_XDIR_,1,solver,mpi);

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  gpuNavierStokes3DParabolicFunction_par_kernel<<<nblocks_par, GPU_THREADS_PER_BLOCK>>>(
      solver->npoints_local, size, ghosts, _XDIR_, dim, dxinv, FDeriv, par
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  memory_accessed = (solver->npoints_local + 2*solver->npoints_local*_MODEL_NVARS_)*sizeof(double);
  printf("%-50s GPU time (secs) = %.6f bandwidth (GB/s) = %6.2f\n",
          "NavierStokes3DParabolicFunction_par", milliseconds*1e-3,
          (1e-6*(memory_accessed)/milliseconds));
#endif

  /* Along Y */

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  gpuNavierStokes3DParabolicFunction_Ydir_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      size, two_third, inv_gamma_m1, inv_Re, inv_Pr, Q,
      QDerivX, QDerivY, QDerivZ, FViscous
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  memory_accessed = 17*size*sizeof(double);
  printf("%-50s GPU time (secs) = %.6f bandwidth (GB/s) = %6.2f\n",
          "NavierStokes3DParabolicFunction_Ydir", milliseconds*1e-3,
          (1e-6*(memory_accessed)/milliseconds));
#endif

  solver->FirstDerivativePar(FDeriv,FViscous,_YDIR_,1,solver,mpi);
  gpuNavierStokes3DParabolicFunction_par_kernel<<<nblocks_par, GPU_THREADS_PER_BLOCK>>>(
      solver->npoints_local, size, ghosts, _YDIR_, dim, dxinv, FDeriv, par
  );
  cudaDeviceSynchronize();

  /* Along Z */
#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(startEvent, 0) );
#endif

  gpuNavierStokes3DParabolicFunction_Zdir_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      size, two_third, inv_gamma_m1, inv_Re, inv_Pr, Q,
      QDerivX, QDerivY, QDerivZ, FViscous
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventRecord(stopEvent, 0) );
  checkCuda( cudaEventSynchronize(stopEvent) );
  checkCuda( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent) );

  memory_accessed = 17*size*sizeof(double);
  printf("%-50s GPU time (secs) = %.6f bandwidth (GB/s) = %6.2f\n",
          "NavierStokes3DParabolicFunction_Zdir", milliseconds*1e-3,
          (1e-6*(memory_accessed)/milliseconds));
#endif

  solver->FirstDerivativePar(FDeriv,FViscous,_ZDIR_,1,solver,mpi);
  gpuNavierStokes3DParabolicFunction_par_kernel<<<nblocks_par, GPU_THREADS_PER_BLOCK>>>(
      solver->npoints_local, size, ghosts, _ZDIR_, dim, dxinv, FDeriv, par
  );
  cudaDeviceSynchronize();

#if defined(GPU_STAT)
  checkCuda( cudaEventDestroy(startEvent) );
  checkCuda( cudaEventDestroy(stopEvent) );
#endif

  if (solver->flag_ib) gpuArrayBlockMultiply(par,solver->gpu_iblank,size,_MODEL_NVARS_);

  return (0);
}

#endif
