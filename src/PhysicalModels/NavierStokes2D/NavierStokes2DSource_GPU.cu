/*! @file NavierStokes2DSource_GPU.cu
    @author Youngdae Kim
    @brief Compute the gravitational source term for the 2D Navier Stokes system
*/
#include <stdlib.h>
#include <basic_gpu.h>
#include <arrayfunctions.h>
#include <physicalmodels/navierstokes2d.h>
#include <mpivars.h>
#include <hypar.h>

#ifdef CUDA_VAR_ORDERDING_AOS

static int gpuNavierStokes2DSourceFunction (double*,const double*,const double*,void*,void*,double,int);
static int gpuNavierStokes2DSourceUpwind   (double*,const double*,const double*,const double*,int,void*,double);

__global__
void NavierStokes2DSource_Xdir_kernel(
  int ngrid_points,
  int ghosts,
  int ndims,
  double gamma,
  double RT,
  const int *dim,
  const double *dxinv,
  const double *grav_field_f,
  const double *SourceI,
  const double *u,
  double *source
)
{
  int tx = threadIdx.x + (blockDim.x * blockIdx.x);
  if (tx < ngrid_points) {
    int v,p,p1,p2;
    int index[GPU_MAX_NDIMS],index1[GPU_MAX_NDIMS],index2[GPU_MAX_NDIMS],dim_interface[GPU_MAX_NDIMS];

    _ArrayIndexnD_(ndims,tx,dim,index,0);
    _ArrayCopy1D_(dim,dim_interface,ndims); dim_interface[_XDIR_]++;
    _ArrayCopy1D_(index,index1,ndims);
    _ArrayCopy1D_(index,index2,ndims); index2[_XDIR_]++;
    _ArrayIndex1D_(ndims,dim          ,index ,ghosts,p );
    _ArrayIndex1D_(ndims,dim_interface,index1,0     ,p1);
    _ArrayIndex1D_(ndims,dim_interface,index2,0     ,p2);

    double dx_inverse; _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,dxinv,dx_inverse);
    double rho, uvel, vvel, e, P; _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*p),rho,uvel,vvel,e,P,gamma);
    double term[_MODEL_NVARS_] = {0.0, rho*RT, 0.0, rho*RT*vvel};
    for (v=0; v<_MODEL_NVARS_; v++) {
      source[_MODEL_NVARS_*p+v] += (  (term[v]*grav_field_f[p])
                                    * (SourceI[_MODEL_NVARS_*p2+v]-SourceI[_MODEL_NVARS_*p1+v])*dx_inverse );
    }
    uvel = P; /* useless statement to avoid compiler warnings */
  }
  return;
}

__global__
void NavierStokes2DSource_Ydir_kernel(
  int ngrid_points,
  int ghosts,
  int ndims,
  double gamma,
  double RT,
  const int *dim,
  const double *dxinv,
  const double *grav_field_f,
  const double *SourceI,
  const double *u,
  double *source
)
{
  int tx = threadIdx.x + (blockDim.x * blockIdx.x);
  if (tx < ngrid_points) {
    int v,p,p1,p2;
    int index[GPU_MAX_NDIMS],index1[GPU_MAX_NDIMS],index2[GPU_MAX_NDIMS],dim_interface[GPU_MAX_NDIMS];

    _ArrayIndexnD_(ndims,tx,dim,index,0);
    _ArrayCopy1D_(dim,dim_interface,ndims); dim_interface[_YDIR_]++;
    _ArrayCopy1D_(index,index1,ndims);
    _ArrayCopy1D_(index,index2,ndims); index2[_YDIR_]++;
    _ArrayIndex1D_(ndims,dim          ,index ,ghosts,p );
    _ArrayIndex1D_(ndims,dim_interface,index1,0     ,p1);
    _ArrayIndex1D_(ndims,dim_interface,index2,0     ,p2);

    double dy_inverse; _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,dxinv,dy_inverse);
    double rho, uvel, vvel, e, P; _NavierStokes2DGetFlowVar_((u+_MODEL_NVARS_*p),rho,uvel,vvel,e,P,gamma);
    double term[_MODEL_NVARS_] = {0.0, 0.0, rho*RT, rho*RT*vvel};
    for (v=0; v<_MODEL_NVARS_; v++) {
      source[_MODEL_NVARS_*p+v] += (  (term[v]*grav_field_f[p])
                                    * (SourceI[_MODEL_NVARS_*p2+v]-SourceI[_MODEL_NVARS_*p1+v])*dy_inverse );
    }
    uvel = P; /* useless statement to avoid compiler warnings */
  }
  return;
}

__global__
void NavierStokes2DSourceUpwind_kernel(
    int ngrid_points,
    const double *fL,
    const double *fR,
    double *fI
)
{
    int p = threadIdx.x + (blockDim.x * blockIdx.x);
    if (p < ngrid_points) {
        for (int k = 0; k < _MODEL_NVARS_; k++)
            (fI+_MODEL_NVARS_*p)[k] = 0.5 * ((fL+_MODEL_NVARS_*p)[k] + (fR+_MODEL_NVARS_*p)[k]);
    }
    return;
}

__global__
void NavierStokes2DSourceFunction_kernel(
    int ngrid_points,
    int dir,
    const double *grav_field_g,
    double *f
)
{
    int p = threadIdx.x + (blockDim.x * blockIdx.x);
    if (p < ngrid_points) {
        (f+_MODEL_NVARS_*p)[0] = 0.0;
        (f+_MODEL_NVARS_*p)[1] = grav_field_g[p] * (dir == _XDIR_);
        (f+_MODEL_NVARS_*p)[2] = grav_field_g[p] * (dir == _YDIR_);
        (f+_MODEL_NVARS_*p)[3] = grav_field_g[p];
    }
    return;
}

/*! Computes the gravitational source term using a well-balanced formulation: the source term
    is rewritten as follows:
    \f{equation}{
      \left[\begin{array}{c} 0 \\ -\rho {\bf g}\cdot{\bf \hat{i}} \\ -\rho {\bf g}\cdot{\bf \hat{j}}  \\ -\rho u {\bf g}\cdot{\bf \hat{i}} - \rho v {\bf g}\cdot{\bf \hat{j}} \end{array}\right]
      =
      \left[\begin{array}{cccc} 0 & p_0 \varrho^{-1} & 0 &  p_0 u \varrho^{-1} \end{array}\right] \cdot \frac{\partial}{\partial x}\left[\begin{array}{c} 0 \\ \varphi \\ 0 \\ \varphi \end{array}\right]
      +
      \left[\begin{array}{cccc} 0 & 0 & p_0 \varrho^{-1} &  p_0 v \varrho^{-1} \end{array}\right] \cdot \frac{\partial}{\partial y}\left[\begin{array}{c} 0 \\ 0 \\ \varphi \\  \varphi \end{array}\right]
    \f}
    where \f$\varphi = \varphi\left(x,y\right)\f$ and \f$\varrho = \varrho\left(x,y\right)\f$ are computed in
    NavierStokes2DGravityField(). The derivatives are computed in an exactly identical manner as the hyperbolic
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
extern "C" int gpuNavierStokes2DSource(
    double  *source,  /*!< Array to hold the computed source */
    double  *gpu_u,       /*!< Solution vector array */
    void    *s,       /*!< Solver object of type #HyPar */
    void    *m,       /*!< MPI object of type #MPIVariables */
    double  t         /*!< Current simulation time */
)
{
  HyPar           *solver = (HyPar* )         s;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;

  if ((param->grav_x == 0.0) && (param->grav_y == 0.0))
    return(0); /* no gravitational forces */

  double  *SourceI = solver->fluxI; /* interace source term       */
  double  *SourceC = solver->fluxC; /* cell-centered source term  */
  double  *SourceL = solver->fL;
  double  *SourceR = solver->fR;

  int     ndims      = solver->ndims;
  int     ghosts     = solver->ghosts;
  int     *dim       = solver->dim_local;
  double  *gpu_x      = solver->gpu_x;
  double  *gpu_dxinv  = solver->gpu_dxinv;
  double  RT         = param->p0 / param->rho0;

  /* Along X-direction */
  if (param->grav_x != 0.0) {
    printf("ALONG X-DIRECTION not tested yet\n");
    exit(0);

    /* calculate the split source function exp(-phi/RT) */
    gpuNavierStokes2DSourceFunction(SourceC,gpu_u,gpu_x,solver,mpi,t,_XDIR_);
    /* calculate the left and right interface source terms */
    solver->InterpolateInterfacesHyp(SourceL,SourceC,gpu_u,gpu_x, 1,_XDIR_,solver,mpi,0);
    solver->InterpolateInterfacesHyp(SourceR,SourceC,gpu_u,gpu_x,-1,_XDIR_,solver,mpi,0);
    /* calculate the final interface source term */
    gpuNavierStokes2DSourceUpwind(SourceI,SourceL,SourceR,gpu_u,_XDIR_,solver,t);
    /* calculate the final cell-centered source term */

    int ngrid_points = 1; for (int i = 0; i < ndims; i++) ngrid_points *= dim[i];
    int nblocks = (ngrid_points - 1) / GPU_THREADS_PER_BLOCK + 1;
    NavierStokes2DSource_Xdir_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      ngrid_points, ghosts, ndims, param->gamma, RT,
      solver->gpu_dim_local, gpu_dxinv, param->gpu_grav_field_f, SourceI, gpu_u, source
    );
  }

  /* Along Y-direction */
  if (param->grav_y != 0.0) {
    printf("ALONG Y-DIRECTION\n");
    /* set interface dimensions */
    /*_ArrayCopy1D_(dim,dim_interface,ndims); dim_interface[_YDIR_]++;*/
    /* calculate the split source function exp(-phi/RT) */
    gpuNavierStokes2DSourceFunction(SourceC,gpu_u,gpu_x,solver,mpi,t,_YDIR_);
    /* calculate the left and right interface source terms */
    solver->InterpolateInterfacesHyp(SourceL,SourceC,gpu_u,gpu_x, 1,_YDIR_,solver,mpi,0);
    solver->InterpolateInterfacesHyp(SourceR,SourceC,gpu_u,gpu_x,-1,_YDIR_,solver,mpi,0);
    /* calculate the final interface source term */
    gpuNavierStokes2DSourceUpwind(SourceI,SourceL,SourceR,gpu_u,_YDIR_,solver,t);

    int ngrid_points = 1; for (int i = 0; i < ndims; i++) ngrid_points *= dim[i];
    int nblocks = (ngrid_points - 1) / GPU_THREADS_PER_BLOCK + 1;
    NavierStokes2DSource_Ydir_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
      ngrid_points, ghosts, ndims, param->gamma, RT,
      solver->gpu_dim_local, gpu_dxinv, param->gpu_grav_field_f, SourceI, gpu_u, source
    );
  }
  return(0);
}

/*! Compute the source function in the well-balanced treatment of the source terms. The source
    function is:
    \f{equation}{
      dir = x \rightarrow \left[\begin{array}{c}0 \\ \varphi\left(x,y\right) \\ 0 \\ \varphi\left(x,y\right) \end{array}\right],
      \ dir = y \rightarrow \left[\begin{array}{c}0 \\ 0 \\ \varphi\left(x,y\right) \\ \varphi\left(x,y\right) \end{array}\right]
    \f}
    where \f$\varphi\f$ (#NavierStokes2D::grav_field_g) is computed in NavierStokes2D::GravityField().
    \n\n
    References:
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, Submitted
*/
int gpuNavierStokes2DSourceFunction(
    double  *gpu_f, /*!< Array to hold the computed source function */
    const double  *gpu_u, /*!< Solution vector array */
    const double  *gpu_x, /*!< Array of spatial coordinates (grid) */
    void    *s, /*!< Solver object of type #HyPar */
    void    *m, /*!< MPI object of type #MPIVariables */
    double  t,  /*!< Current simulation time */
    int     dir /*!< Spatial dimension (x or y) */
)
{
  HyPar           *solver = (HyPar* )         s;
  NavierStokes2D  *param  = (NavierStokes2D*) solver->physics;

  int     ghosts  = solver->ghosts;
  int     *dim    = solver->dim_local;
  int     ndims   = solver->ndims;

  double cpu_time = 0.0;
  clock_t cpu_start, cpu_end;

  int ngrid_points = 1; for (int i = 0; i < ndims; i++) ngrid_points *= (dim[i] + 2*ghosts);
  int nblocks = (ngrid_points - 1) / GPU_THREADS_PER_BLOCK + 1;

  cpu_start = clock();
  NavierStokes2DSourceFunction_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    ngrid_points, dir, param->gpu_grav_field_g, gpu_f);
  cudaDeviceSynchronize();
  cpu_end = clock();
  cpu_time += (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;

  return(0);
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
int gpuNavierStokes2DSourceUpwind(
  double  *gpu_fI,  /*!< Array to hold the computed "upwind" interface source function */
  const double  *gpu_fL,  /*!< Interface source function value computed using left-biased interpolation */
  const double  *gpu_fR,  /*!< Interface source function value computed using right-biased interpolation */
  const double  *gpu_u,   /*!< Solution vector array */
  int     dir,  /*!< Spatial dimension (x or y) */
  void    *s,   /*!< Solver object of type #HyPar */
  double  t     /*!< Current simulation time */
)
{
  HyPar           *solver = (HyPar*) s;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int *dim  = solver->dim_local;
  int bounds_inter[ndims];

  double cpu_time = 0.0;
  clock_t cpu_start, cpu_end;

  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  int ngrid_points; _ArrayProduct1D_(bounds_inter,ndims,ngrid_points);
  int nblocks = (ngrid_points - 1) / GPU_THREADS_PER_BLOCK + 1;

  cpu_start = clock();
  NavierStokes2DSourceUpwind_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
    ngrid_points, gpu_fL, gpu_fR, gpu_fI);
  cudaDeviceSynchronize();
  cpu_end = clock();
  cpu_time += (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;

  return(0);
}

#endif
