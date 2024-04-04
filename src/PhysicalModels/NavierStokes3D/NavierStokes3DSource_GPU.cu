/*! @file NavierStokes3DSource_GPU.cu
    @author Youngdae Kim
    @brief Compute the gravitational source term for the 3D Navier Stokes system
*/
#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>
#include <physicalmodels/navierstokes3d.h>
#include <mpivars.h>
#include <hypar.h>

static int gpuNavierStokes3DSourceFunction (double*,double*,double*,void*,void*,double,int);
static int gpuNavierStokes3DSourceUpwind   (double*,double*,double*,double*,int,void*,double);

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel function */
__global__
void gpuNavierStokes3DSource_grav_kernel(
  int npoints_grid,
  int ghosts,
  int dir,
  double gamma,
  double RT,
  const int * __restrict__ dim,
  const double * __restrict__ dxinv,
  const double * __restrict__ grav_field_f,
  const double * __restrict__ SourceI,
  const double * __restrict__ u,
  double * __restrict__ source
)
{
  int tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < npoints_grid) {
    int v,p,p1,p2;
    int index[GPU_MAX_NDIMS],index1[GPU_MAX_NDIMS],index2[GPU_MAX_NDIMS],dim_interface[GPU_MAX_NDIMS];

    _ArrayIndexnD_(_MODEL_NDIMS_,tid,dim,index,0);
    _ArrayCopy1D_(dim,dim_interface,_MODEL_NDIMS_); dim_interface[dir]++;
    _ArrayCopy1D_(index,index1,_MODEL_NDIMS_);
    _ArrayCopy1D_(index,index2,_MODEL_NDIMS_); index2[dir]++;
    _ArrayIndex1D_(_MODEL_NDIMS_,dim          ,index ,ghosts,p );
    _ArrayIndex1D_(_MODEL_NDIMS_,dim_interface,index1,0     ,p1);
    _ArrayIndex1D_(_MODEL_NDIMS_,dim_interface,index2,0     ,p2);

    double dx_inverse; _GetCoordinate_(dir,index[dir],dim,ghosts,dxinv,dx_inverse);
    double rho, vel[_MODEL_NDIMS_], e, P;
    _NavierStokes3DGetFlowVar_((u+_MODEL_NVARS_*p),_NavierStokes3D_stride_,rho,vel[0],vel[1],vel[2],e,P,gamma);
    double term[_MODEL_NVARS_] = {0.0, rho*RT*(dir==_XDIR_), rho*RT*(dir==_YDIR_), rho*RT*(dir==_ZDIR_), rho*RT*vel[dir]};
    for (v=0; v<_MODEL_NVARS_; v++) {
      source[_MODEL_NVARS_*p+v] += (  (term[v]*grav_field_f[p])
                                    * (SourceI[_MODEL_NVARS_*p2+v]-SourceI[_MODEL_NVARS_*p1+v])*dx_inverse );
    }
  }
  return;
}

/*! Kernel function for gpuNavierStokes3DSourceFunction() */
__global__
void gpuNavierStokes3DSourceFunction_kernel(
  int npoints_grid,
  int dir,
  const double * __restrict__ grav_field_g,
  double * __restrict__ f
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    (f+_MODEL_NVARS_*p)[0] = 0.0;
    (f+_MODEL_NVARS_*p)[1] = grav_field_g[p] * (dir == _XDIR_);
    (f+_MODEL_NVARS_*p)[2] = grav_field_g[p] * (dir == _YDIR_);
    (f+_MODEL_NVARS_*p)[3] = grav_field_g[p] * (dir == _ZDIR_);
    (f+_MODEL_NVARS_*p)[4] = grav_field_g[p];
  }
  return;
}

/*! Kernel function for gpuNavierStokes3DSourceUpwind() */
__global__
void gpuNavierStokes3DSourceUpwind_kernel(
  int npoints_grid,
  const double * __restrict__ fL,
  const double * __restrict__ fR,
  double * __restrict__ fI
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
      for (int k = 0; k < _MODEL_NVARS_; k++)
        (fI+_MODEL_NVARS_*p)[k] = 0.5 * ((fL+_MODEL_NVARS_*p)[k] + (fR+_MODEL_NVARS_*p)[k]);
  }
  return;
}

/*! Compute the source function in the well-balanced treatment of the source terms. The source
    function is:
    \f{equation}{
      dir = x \rightarrow \left[\begin{array}{c}0 \\ \varphi\left(x,y,z\right) \\ 0 \\ 0 \\ \varphi\left(x,y,z\right) \end{array}\right],
      \ dir = y \rightarrow \left[\begin{array}{c}0 \\ 0 \\ \varphi\left(x,y,z\right) \\ 0 \\ \varphi\left(x,y,z\right) \end{array}\right],
      \ dir = z \rightarrow \left[\begin{array}{c}0 \\ 0 \\ 0 \\ \varphi\left(x,y,z\right) \\ \varphi\left(x,y,z\right) \end{array}\right]
    \f}
    where \f$\varphi\f$ (#NavierStokes3D::grav_field_g) is computed in NavierStokes3D::GravityField().
    \n\n
    References:
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, AIAA Journal, http://dx.doi.org/10.2514/1.J054580.
*/
int gpuNavierStokes3DSourceFunction(
  double  *f, /*!< Array to hold the computed source function */
  double  *u, /*!< Solution vector array */
  double  *x, /*!< Array of spatial coordinates (grid) */
  void    *s, /*!< Solver object of type #HyPar */
  void    *m, /*!< MPI object of type #MPIVariables */
  double  t,  /*!< Current simulation time */
  int     dir /*!< Spatial dimension (x, y, or z) */
)
{
  HyPar           *solver = (HyPar* )         s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;
  int             nblocks = (solver->npoints_local_wghosts-1)/GPU_THREADS_PER_BLOCK + 1;

  gpuNavierStokes3DSourceFunction_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    solver->npoints_local_wghosts, dir, param->gpu_grav_field_g, f
  );
  cudaDeviceSynchronize();

  return 0;
}

/*! Compute the "upwind" source function value at the interface: the upwinding is just the
    arithmetic average of the left and right biased interpolated values of the source function
    at the interface.
    \n\n
    References:
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, Submitted
*/
int gpuNavierStokes3DSourceUpwind(
  double  *fI,  /*!< Array to hold the computed "upwind" interface source function */
  double  *fL,  /*!< Interface source function value computed using left-biased interpolation */
  double  *fR,  /*!< Interface source function value computed using right-biased interpolation */
  double  *u,   /*!< Solution vector array */
  int     dir,  /*!< Spatial dimension (x,y, or z) */
  void    *s,   /*!< Solver object of type #HyPar */
  double  t     /*!< Current simulation time */
)
{
  HyPar *solver = (HyPar*) s;
  int   *dim    = solver->dim_local;
  _DECLARE_IERR_;

  int bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;
  int npoints_grid; _ArrayProduct1D_(bounds_inter,_MODEL_NDIMS_,npoints_grid);
  int nblocks = (npoints_grid-1)/GPU_THREADS_PER_BLOCK + 1;

  gpuNavierStokes3DSourceUpwind_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    npoints_grid, fL, fR, fI
  );
  cudaDeviceSynchronize();

  return 0;
}

/*! Computes the gravitational source term using a well-balanced formulation: the source term
    is rewritten as follows:
    \f{equation}{
      \left[\begin{array}{c} 0 \\ -\rho {\bf g}\cdot{\bf \hat{i}} \\ -\rho {\bf g}\cdot{\bf \hat{j}} \\ -\rho {\bf g}\cdot{\bf \hat{k}}  \\ -\rho u {\bf g}\cdot{\bf \hat{i}} - \rho v {\bf g}\cdot{\bf \hat{j}} - \rho w {\bf g}\cdot{\bf \hat{k}} \end{array}\right]
      =
      \left[\begin{array}{ccccc} 0 & p_0 \varrho^{-1} & 0 & 0 &  p_0 u \varrho^{-1} \end{array}\right] \cdot \frac{\partial}{\partial x}\left[\begin{array}{c} 0 \\ \varphi \\ 0 \\ 0 \\ \varphi \end{array}\right]
      +
      \left[\begin{array}{ccccc} 0 & 0 & p_0 \varrho^{-1} & 0 &  p_0 v \varrho^{-1} \end{array}\right] \cdot \frac{\partial}{\partial y}\left[\begin{array}{c} 0 \\ 0 \\ \varphi \\ 0 \\ \varphi \end{array}\right]
      +
      \left[\begin{array}{ccccc} 0 & 0 & 0 &  p_0 \varrho^{-1} & p_0 w \varrho^{-1} \end{array}\right] \cdot \frac{\partial}{\partial y}\left[\begin{array}{c} 0 \\ 0 \\ 0 \\ \varphi \\ \varphi \end{array}\right]
    \f}
    where \f$\varphi = \varphi\left(x,y\right)\f$ and \f$\varrho = \varrho\left(x,y\right)\f$ are computed in
    NavierStokes3DGravityField(). The derivatives are computed in an exactly identical manner as the hyperbolic
    flux (i.e., using a conservative finite-difference formulation) (see HyperbolicFunction()).
    \n\n
    References:
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, AIAA Journal, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580.
*/
extern "C" int gpuNavierStokes3DSource(
  double * __restrict__ source, /*!< Array to hold the computed source */
  double * __restrict__ u,      /*!< Solution vector array */
  void   * __restrict__ s,      /*!< Solver object of type #HyPar */
  void   * __restrict__ m,      /*!< MPI object of type #MPIVariables */
  double t                      /*!< Current simulation time */
)
{
  HyPar           *solver = (HyPar* )         s;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;

  if ((param->grav_x == 0.0) && (param->grav_y == 0.0) && (param->grav_z == 0.0))
    return(0); /* no gravitational forces */

  int     dir;
  double  *SourceI = solver->fluxI; /* interace source term       */
  double  *SourceC = solver->fluxC; /* cell-centered source term  */
  double  *SourceL = solver->fL;
  double  *SourceR = solver->fR;

  int     ghosts  = solver->ghosts;
  int     *dim    = solver->gpu_dim_local;
  double  *x      = solver->gpu_x;
  double  *dxinv  = solver->gpu_dxinv;
  double  RT      =  param->p0 / param->rho0;
  static double grav[_MODEL_NDIMS_];

  grav[_XDIR_] = param->grav_x;
  grav[_YDIR_] = param->grav_y;
  grav[_ZDIR_] = param->grav_z;

  int nblocks = (solver->npoints_local-1)/GPU_THREADS_PER_BLOCK + 1;
  for (dir = 0; dir < _MODEL_NDIMS_; dir++) {
    if (grav[dir] != 0.0) {
      /* calculate the split source function exp(-phi/RT) */
      gpuNavierStokes3DSourceFunction(SourceC,u,x,solver,mpi,t,dir);
      /* calculate the left and right interface source terms */
      solver->InterpolateInterfacesHyp(SourceL,SourceC,u,x, 1,dir,solver,mpi,0);
      solver->InterpolateInterfacesHyp(SourceR,SourceC,u,x,-1,dir,solver,mpi,0);
      /* calculate the final interface source term */
      gpuNavierStokes3DSourceUpwind(SourceI,SourceL,SourceR,u,dir,solver,t);
      /* calculate the final cell-centered source term */
      gpuNavierStokes3DSource_grav_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        solver->npoints_local, ghosts, dir, param->gamma, RT,
        dim, dxinv, param->gpu_grav_field_f, SourceI, u, source
      );
    }
  }

  return 0;
}

#else

/*! Kernel function */
__global__
void gpuNavierStokes3DSource_grav_kernel(
  int npoints_grid,
  int npoints_local_wghosts,
  int npoints_fluxI,
  int ghosts,
  int dir,
  double gamma,
  double RT,
  const int * __restrict__ dim,
  const double * __restrict__ dxinv,
  const double * __restrict__ grav_field_f,
  const double * __restrict__ SourceI,
  const double * __restrict__ u,
  double * __restrict__ source
)
{
  int tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < npoints_grid) {
    int v,p,p1,p2;
    int index[_MODEL_NDIMS_],index1[_MODEL_NDIMS_],index2[_MODEL_NDIMS_],dim_interface[_MODEL_NDIMS_];

    _ArrayIndexnD_(_MODEL_NDIMS_,tid,dim,index,0);
    _ArrayCopy1D_(dim,dim_interface,_MODEL_NDIMS_); dim_interface[dir]++;
    _ArrayCopy1D_(index,index1,_MODEL_NDIMS_);
    _ArrayCopy1D_(index,index2,_MODEL_NDIMS_); index2[dir]++;
    _ArrayIndex1D_(_MODEL_NDIMS_,dim          ,index ,ghosts,p );
    _ArrayIndex1D_(_MODEL_NDIMS_,dim_interface,index1,0     ,p1);
    _ArrayIndex1D_(_MODEL_NDIMS_,dim_interface,index2,0     ,p2);

    double dx_inverse; _GetCoordinate_(dir,index[dir],dim,ghosts,dxinv,dx_inverse);
    double rho, vel[_MODEL_NDIMS_], e, P;
    _NavierStokes3DGetFlowVar_((u+p),npoints_local_wghosts,rho,vel[0],vel[1],vel[2],e,P,gamma);
    double term[_MODEL_NVARS_] = {0.0, rho*RT*(dir==_XDIR_), rho*RT*(dir==_YDIR_), rho*RT*(dir==_ZDIR_), rho*RT*vel[dir]};
    for (v=0; v<_MODEL_NVARS_; v++) {
      source[p+v*npoints_local_wghosts] += (  (term[v]*grav_field_f[p])
                                            * (SourceI[p2+v*npoints_fluxI]-SourceI[p1+v*npoints_fluxI])*dx_inverse );
    }
    vel[0] = P;
  }
  return;
}

/*! Kernel function for gpuNavierStokes3DSourceFunction() */
__global__
void gpuNavierStokes3DSourceFunction_kernel(
  int npoints_grid,
  int ghosts,
  int dir,
  const int * __restrict__ dim,
  const double * __restrict__ grav_field_g,
  double * __restrict__ f
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    (f+p)[0*npoints_grid] = 0.0;
    (f+p)[1*npoints_grid] = grav_field_g[p] * (dir == _XDIR_);
    (f+p)[2*npoints_grid] = grav_field_g[p] * (dir == _YDIR_);
    (f+p)[3*npoints_grid] = grav_field_g[p] * (dir == _ZDIR_);
    (f+p)[4*npoints_grid] = grav_field_g[p];
  }
  return;
}

/*! Kernel function for gpuNavierStokes3DSourceUpwind() */
__global__
void gpuNavierStokes3DSourceUpwind_kernel(
  int npoints_grid,
  const double * __restrict__ fL,
  const double * __restrict__ fR,
  double * __restrict__ fI
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
      for (int k = 0; k < _MODEL_NVARS_; k++)
        (fI+p)[k*npoints_grid] = 0.5 * ((fL+p)[k*npoints_grid] + (fR+p)[k*npoints_grid]);
  }
  return;
}

/*! Compute the source function in the well-balanced treatment of the source terms. The source
    function is:
    \f{equation}{
      dir = x \rightarrow \left[\begin{array}{c}0 \\ \varphi\left(x,y,z\right) \\ 0 \\ 0 \\ \varphi\left(x,y,z\right) \end{array}\right],
      \ dir = y \rightarrow \left[\begin{array}{c}0 \\ 0 \\ \varphi\left(x,y,z\right) \\ 0 \\ \varphi\left(x,y,z\right) \end{array}\right],
      \ dir = z \rightarrow \left[\begin{array}{c}0 \\ 0 \\ 0 \\ \varphi\left(x,y,z\right) \\ \varphi\left(x,y,z\right) \end{array}\right]
    \f}
    where \f$\varphi\f$ (#NavierStokes3D::grav_field_g) is computed in NavierStokes3D::GravityField().
    \n\n
    References:
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, AIAA Journal, http://dx.doi.org/10.2514/1.J054580.
*/
int gpuNavierStokes3DSourceFunction(
  double  *f, /*!< Array to hold the computed source function */
  double  *u, /*!< Solution vector array */
  double  *x, /*!< Array of spatial coordinates (grid) */
  void    *s, /*!< Solver object of type #HyPar */
  void    *m, /*!< MPI object of type #MPIVariables */
  double  t,  /*!< Current simulation time */
  int     dir /*!< Spatial dimension (x, y, or z) */
)
{
  HyPar           *solver = (HyPar* )         s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;
  int             nblocks = (solver->npoints_local_wghosts-1)/GPU_THREADS_PER_BLOCK + 1;

  gpuNavierStokes3DSourceFunction_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    solver->npoints_local_wghosts, solver->ghosts, dir, solver->gpu_dim_local, param->gpu_grav_field_g, f
  );
  cudaDeviceSynchronize();

  return 0;
}

/*! Compute the "upwind" source function value at the interface: the upwinding is just the
    arithmetic average of the left and right biased interpolated values of the source function
    at the interface.
    \n\n
    References:
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, Submitted
*/
int gpuNavierStokes3DSourceUpwind(
  double  *fI,  /*!< Array to hold the computed "upwind" interface source function */
  double  *fL,  /*!< Interface source function value computed using left-biased interpolation */
  double  *fR,  /*!< Interface source function value computed using right-biased interpolation */
  double  *u,   /*!< Solution vector array */
  int     dir,  /*!< Spatial dimension (x,y, or z) */
  void    *s,   /*!< Solver object of type #HyPar */
  double  t     /*!< Current simulation time */
)
{
  HyPar *solver = (HyPar*) s;
  int   *dim    = solver->dim_local;

  int bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;
  int npoints_grid; _ArrayProduct1D_(bounds_inter,_MODEL_NDIMS_,npoints_grid);
  int nblocks = (npoints_grid-1)/GPU_THREADS_PER_BLOCK + 1;

  gpuNavierStokes3DSourceUpwind_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    npoints_grid, fL, fR, fI
  );
  cudaDeviceSynchronize();

  return 0;
}

/*! Computes the gravitational source term using a well-balanced formulation: the source term
    is rewritten as follows:
    \f{equation}{
      \left[\begin{array}{c} 0 \\ -\rho {\bf g}\cdot{\bf \hat{i}} \\ -\rho {\bf g}\cdot{\bf \hat{j}} \\ -\rho {\bf g}\cdot{\bf \hat{k}}  \\ -\rho u {\bf g}\cdot{\bf \hat{i}} - \rho v {\bf g}\cdot{\bf \hat{j}} - \rho w {\bf g}\cdot{\bf \hat{k}} \end{array}\right]
      =
      \left[\begin{array}{ccccc} 0 & p_0 \varrho^{-1} & 0 & 0 &  p_0 u \varrho^{-1} \end{array}\right] \cdot \frac{\partial}{\partial x}\left[\begin{array}{c} 0 \\ \varphi \\ 0 \\ 0 \\ \varphi \end{array}\right]
      +
      \left[\begin{array}{ccccc} 0 & 0 & p_0 \varrho^{-1} & 0 &  p_0 v \varrho^{-1} \end{array}\right] \cdot \frac{\partial}{\partial y}\left[\begin{array}{c} 0 \\ 0 \\ \varphi \\ 0 \\ \varphi \end{array}\right]
      +
      \left[\begin{array}{ccccc} 0 & 0 & 0 &  p_0 \varrho^{-1} & p_0 w \varrho^{-1} \end{array}\right] \cdot \frac{\partial}{\partial y}\left[\begin{array}{c} 0 \\ 0 \\ 0 \\ \varphi \\ \varphi \end{array}\right]
    \f}
    where \f$\varphi = \varphi\left(x,y\right)\f$ and \f$\varrho = \varrho\left(x,y\right)\f$ are computed in
    NavierStokes3DGravityField(). The derivatives are computed in an exactly identical manner as the hyperbolic
    flux (i.e., using a conservative finite-difference formulation) (see HyperbolicFunction()).
    \n\n
    References:
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, AIAA Journal, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580.
*/
extern "C" int gpuNavierStokes3DSource(
  double * __restrict__ source,
  double * __restrict__ u,
  void   * __restrict__ s,
  void   * __restrict__ m,
  double t
)
{
  HyPar           *solver = (HyPar* )         s;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;

  if ((param->grav_x == 0.0) && (param->grav_y == 0.0) && (param->grav_z == 0.0))
    return(0); /* no gravitational forces */

  int     dir;
  double  *SourceI = solver->fluxI; /* interace source term       */
  double  *SourceC = solver->fluxC; /* cell-centered source term  */
  double  *SourceL = solver->fL;
  double  *SourceR = solver->fR;

  int     ghosts  = solver->ghosts;
  int     *dim    = solver->gpu_dim_local;
  double  *x      = solver->gpu_x;
  double  *dxinv  = solver->gpu_dxinv;
  double  RT      = param->p0 / param->rho0;
  static double grav[_MODEL_NDIMS_];

  grav[_XDIR_] = param->grav_x;
  grav[_YDIR_] = param->grav_y;
  grav[_ZDIR_] = param->grav_z;

  int nblocks = (solver->npoints_local-1)/GPU_THREADS_PER_BLOCK + 1;
  for (dir = 0; dir < _MODEL_NDIMS_; dir++) {
    if (grav[dir] != 0.0) {
      int bounds_inter[_MODEL_NDIMS_];
      _ArrayCopy1D_(solver->dim_local,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;
      int npoints_fluxI; _ArrayProduct1D_(bounds_inter,_MODEL_NDIMS_,npoints_fluxI);

      /* calculate the split source function exp(-phi/RT) */
      gpuNavierStokes3DSourceFunction(SourceC,u,x,solver,mpi,t,dir);
      /* calculate the left and right interface source terms */
      solver->InterpolateInterfacesHyp(SourceL,SourceC,u,x, 1,dir,solver,mpi,0);
      solver->InterpolateInterfacesHyp(SourceR,SourceC,u,x,-1,dir,solver,mpi,0);
      /* calculate the final interface source term */
      gpuNavierStokes3DSourceUpwind(SourceI,SourceL,SourceR,u,dir,solver,t);
      /* calculate the final cell-centered source term */
      gpuNavierStokes3DSource_grav_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        solver->npoints_local, solver->npoints_local_wghosts, npoints_fluxI,
        ghosts, dir, param->gamma, RT,
        dim, dxinv, param->gpu_grav_field_f, SourceI, u, source
      );
      cudaDeviceSynchronize();
    }
  }

  return 0;
}

#endif

