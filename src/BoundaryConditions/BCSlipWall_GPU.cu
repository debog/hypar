/*! @file BCSlipWall_GPU.cu
    @author Youngdae Kim
    @brief GPU implmentation of slip-wall boundary conditions
*/
#include <stdlib.h>
#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>
#include <boundaryconditions.h>

#include <physicalmodels/navierstokes2d.h>
#include <physicalmodels/navierstokes3d.h>

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel for 2D implementation of gpuBCSlipWallU() */
__global__
void BCSlipWallU_dim2_kernel(
  int grid_size,
  int dim,
  int face,
  int ghosts,
  int ndims,
  int nvars,
  double gamma,
  const int *size,
  const int *ie,
  const int *is,
  const double *FlowVelocity,
  double *phi
)
{
  int index = threadIdx.x + (blockDim.x * blockIdx.x);
  if (index < grid_size) {
    const int max_ndims = 3;
    int p1, p2;
    int bounds[max_ndims], indexb[max_ndims], indexi[max_ndims];
    double inv_gamma_m1 = 1.0/(gamma-1.0);

    _ArraySubtract1D_(bounds,ie,is,ndims);
    _ArrayIndexnD_(ndims,index,bounds,indexb,0);
    _ArrayCopy1D_(indexb,indexi,ndims);
    _ArrayAdd1D_(indexi,indexi,is,ndims);
    if      (face ==  1) indexi[dim] = ghosts-1-indexb[dim];
    else if (face == -1) indexi[dim] = size[dim]-indexb[dim]-1;
    else return;
    _ArrayIndex1DWO_(ndims,size,indexb,is,ghosts,p1);
    _ArrayIndex1D_(ndims,size,indexi,ghosts,p2);

    /* flow variables in the interior */
    double rho, uvel, vvel, energy, pressure;
    double rho_gpt, uvel_gpt, vvel_gpt, energy_gpt, pressure_gpt;
    _NavierStokes2DGetFlowVar_((phi+nvars*p2),rho,uvel,vvel,energy,pressure,gamma);
    /* set the ghost point values */
    rho_gpt = rho;
    pressure_gpt = pressure;
    if (dim == _XDIR_) {
      uvel_gpt = 2.0*FlowVelocity[_XDIR_] - uvel;
      vvel_gpt = vvel;
    } else if (dim == _YDIR_) {
      uvel_gpt = uvel;
      vvel_gpt = 2.0*FlowVelocity[_YDIR_] - vvel;
    } else {
      uvel_gpt = 0.0;
      vvel_gpt = 0.0;
    }
    energy_gpt = inv_gamma_m1*pressure_gpt
                + 0.5 * rho_gpt * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt);

    phi[nvars*p1+0] = rho_gpt;
    phi[nvars*p1+1] = rho_gpt * uvel_gpt;
    phi[nvars*p1+2] = rho_gpt * vvel_gpt;
    phi[nvars*p1+3] = energy_gpt;
  }
  return;
}

/*! Kernel for 3D implementation of gpuBCSlipWallU() */
__global__
void BCSlipWallU_dim3_kernel(
  int grid_size,
  int dim,
  int face,
  int ghosts,
  int ndims,
  int nvars,
  double gamma,
  const int *size,
  const int *ie,
  const int *is,
  const double *FlowVelocity,
  double *phi
)
{
  int index = threadIdx.x + (blockDim.x * blockIdx.x);
  if (index < grid_size) {
    const int max_ndims = 3;
    int p1, p2;
    int bounds[max_ndims], indexb[max_ndims], indexi[max_ndims];
    double inv_gamma_m1 = 1.0/(gamma-1.0);

    _ArraySubtract1D_(bounds,ie,is,ndims);
    _ArrayIndexnD_(ndims,index,bounds,indexb,0);
    _ArrayCopy1D_(indexb,indexi,ndims);
    _ArrayAdd1D_(indexi,indexi,is,ndims);
    if      (face ==  1) indexi[dim] = ghosts-1-indexb[dim];
    else if (face == -1) indexi[dim] = size[dim]-indexb[dim]-1;
    else return;
    _ArrayIndex1DWO_(ndims,size,indexb,is,ghosts,p1);
    _ArrayIndex1D_(ndims,size,indexi,ghosts,p2);

    /* flow variables in the interior */
    double rho, uvel, vvel, wvel, energy, pressure;
    double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt;
    _NavierStokes3DGetFlowVar_((phi+nvars*p2),_NavierStokes3D_stride_,rho,uvel,vvel,wvel,energy,pressure,gamma);
    /* set the ghost point values */
    rho_gpt = rho;
    pressure_gpt = pressure;
    if (dim == _XDIR_) {
      uvel_gpt = 2.0*FlowVelocity[_XDIR_] - uvel;
      vvel_gpt = vvel;
      wvel_gpt = wvel;
    } else if (dim == _YDIR_) {
      uvel_gpt = uvel;
      vvel_gpt = 2.0*FlowVelocity[_YDIR_] - vvel;
      wvel_gpt = wvel;
    } else if (dim == _ZDIR_) {
      uvel_gpt = uvel;
      vvel_gpt = vvel;
      wvel_gpt = 2.0*FlowVelocity[_ZDIR_] - wvel;
    } else {
      uvel_gpt = 0.0;
      vvel_gpt = 0.0;
      wvel_gpt = 0.0;
    }
    energy_gpt = inv_gamma_m1*pressure_gpt
                + 0.5 * rho_gpt
                * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt + wvel_gpt*wvel_gpt);

    phi[nvars*p1+0] = rho_gpt;
    phi[nvars*p1+1] = rho_gpt * uvel_gpt;
    phi[nvars*p1+2] = rho_gpt * vvel_gpt;
    phi[nvars*p1+3] = rho_gpt * wvel_gpt;
    phi[nvars*p1+4] = energy_gpt;
  }
  return;
}

/*! Applies the slip-wall boundary condition: This is specific to the two and three
    dimensional Euler and Navier-Stokes systems (#Euler2D, #NavierStokes2D, #NavierStokes3D).
    It is used for simulating inviscid walls or symmetric boundaries. The pressure, density,
    and tangential velocity at the ghost points are extrapolated from the interior, while the
    normal velocity at the ghost points is set such that the interpolated value at the boundary
    face is equal to the specified wall velocity.

    \sa BCSlipWallU()
*/
extern "C" int gpuBCSlipWallU(
                 void    *b,     /*!< Boundary object of type #DomainBoundary */
                 void    *m,     /*!< MPI object of type #MPIVariables */
                 int     ndims,  /*!< Number of spatial dimensions */
                 int     nvars,  /*!< Number of variables/DoFs per grid point */
                 int     *gpu_size,  /*!< Integer array with the number of grid points in each spatial dimension */
                 int     ghosts, /*!< Number of ghost points */
                 double  *gpu_phi,   /*!< The solution array on which to apply the boundary condition */
                 double  waqt    /*!< Current solution time */
               )
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;
  int face  = boundary->face;

  if (ndims == 2) {

    if (boundary->on_this_proc) {
      int bounds[ndims];
      _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
      int ngrid_points = 1; for(int i = 0; i < ndims; i++) ngrid_points *= bounds[i];
      int nblocks = (ngrid_points - 1) / GPU_THREADS_PER_BLOCK + 1;

      gpuMemcpy(  boundary->gpu_ie, boundary->ie, ndims*sizeof(int), gpuMemcpyHostToDevice  );
      gpuMemcpy(  boundary->gpu_is, boundary->is, ndims*sizeof(int), gpuMemcpyHostToDevice  );

      BCSlipWallU_dim2_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        ngrid_points, dim, face, ghosts, ndims, nvars, boundary->gamma,
        gpu_size, boundary->gpu_ie, boundary->gpu_is, boundary->gpu_FlowVelocity, gpu_phi);
    }

  } else if (ndims == 3) {

    if (boundary->on_this_proc) {
      int bounds[ndims];
      _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
      int ngrid_points = 1; for(int i = 0; i < ndims; i++) ngrid_points *= bounds[i];
      int nblocks = (ngrid_points - 1) / GPU_THREADS_PER_BLOCK + 1;

      BCSlipWallU_dim3_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        ngrid_points, dim, face, ghosts, ndims, nvars, boundary->gamma,
        gpu_size, boundary->gpu_ie, boundary->gpu_is, boundary->gpu_FlowVelocity, gpu_phi);
    }

  } else {

    fprintf(stderr,"Error in gpuBCSlipWallU(): not implemented for ndims=%d.\n", ndims);
    return 1;

  }

  return 0;
}

#else

/*! Kernel for 3D implementation of gpuBCSlipWallU() */
__global__
void BCSlipWallU_dim3_kernel(
  int npoints_bounds,
  int npoints_local_wghosts,
  int face,
  int ndims,
  int dim,
  int ghosts,
  int nvars,
  double gamma,
  const int * __restrict__ bounds,
  const int * __restrict__ size,
  const int * __restrict__ boundary_is,
  const double * __restrict__ FlowVelocity,
  double * __restrict__ phi
)
{
  int index = threadIdx.x + (blockDim.x * blockIdx.x);
  if (index < npoints_bounds) {
    const int max_ndims = 3;
    int p1, p2;
    int indexb[max_ndims], indexi[max_ndims];
    double inv_gamma_m1 = 1.0/(gamma-1.0);

    _ArrayIndexnD_(ndims,index,bounds,indexb,0);
    _ArrayCopy1D_(indexb,indexi,ndims);
    _ArrayAdd1D_(indexi,indexi,boundary_is,ndims);
    if      (face ==  1) indexi[dim] = ghosts-1-indexb[dim];
    else if (face == -1) indexi[dim] = size[dim]-indexb[dim]-1;
    else return;
    _ArrayIndex1DWO_(ndims,size,indexb,boundary_is,ghosts,p1);
    _ArrayIndex1D_(ndims,size,indexi,ghosts,p2);

    /* flow variables in the interior */
    double rho, uvel, vvel, wvel, energy, pressure;
    double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt;
    _NavierStokes3DGetFlowVar_((phi+p2),npoints_local_wghosts,rho,uvel,vvel,wvel,energy,pressure,gamma);
    /* set the ghost point values */
    rho_gpt = rho;
    pressure_gpt = pressure;
    if (dim == _XDIR_) {
      uvel_gpt = 2.0*FlowVelocity[_XDIR_] - uvel;
      vvel_gpt = vvel;
      wvel_gpt = wvel;
    } else if (dim == _YDIR_) {
      uvel_gpt = uvel;
      vvel_gpt = 2.0*FlowVelocity[_YDIR_] - vvel;
      wvel_gpt = wvel;
    } else if (dim == _ZDIR_) {
      uvel_gpt = uvel;
      vvel_gpt = vvel;
      wvel_gpt = 2.0*FlowVelocity[_ZDIR_] - wvel;
    } else {
      uvel_gpt = 0.0;
      vvel_gpt = 0.0;
      wvel_gpt = 0.0;
    }
    energy_gpt = inv_gamma_m1*pressure_gpt
                + 0.5 * rho_gpt
                * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt + wvel_gpt*wvel_gpt);

    phi[p1+0] = rho_gpt;
    phi[p1+1*npoints_local_wghosts] = rho_gpt * uvel_gpt;
    phi[p1+2*npoints_local_wghosts] = rho_gpt * vvel_gpt;
    phi[p1+3*npoints_local_wghosts] = rho_gpt * wvel_gpt;
    phi[p1+4*npoints_local_wghosts] = energy_gpt;
  }
  return;
}

/*! Applies the slip-wall boundary condition: This is specific to the two and three
    dimensional Euler and Navier-Stokes systems (#Euler2D, #NavierStokes2D, #NavierStokes3D).
    It is used for simulating inviscid walls or symmetric boundaries. The pressure, density,
    and tangential velocity at the ghost points are extrapolated from the interior, while the
    normal velocity at the ghost points is set such that the interpolated value at the boundary
    face is equal to the specified wall velocity.

    \sa BCSlipWallU()
*/
extern "C" int gpuBCSlipWallU(
                 void    * __restrict__ b,     /*!< Boundary object of type #DomainBoundary */
                 void    * __restrict__ m,     /*!< MPI object of type #MPIVariables */
                 int     ndims,  /*!< Number of spatial dimensions */
                 int     nvars,  /*!< Number of variables/DoFs per grid point */
                 int     * __restrict__ size,  /*!< Integer array with the number of grid points in each spatial dimension */
                 int     ghosts, /*!< Number of ghost points */
                 double  * __restrict__ phi,   /*!< The solution array on which to apply the boundary condition */
                 double  waqt    /*!< Current solution time */
               )
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;
  int face  = boundary->face;

  if (ndims == 3) {

    if (boundary->on_this_proc) {
      int nblocks = (boundary->gpu_npoints_bounds - 1) / GPU_THREADS_PER_BLOCK + 1;

      BCSlipWallU_dim3_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        boundary->gpu_npoints_bounds, boundary->gpu_npoints_local_wghosts,
        face, ndims, dim, ghosts, nvars, boundary->gamma,
        boundary->gpu_bounds, size, boundary->gpu_is, boundary->gpu_FlowVelocity, phi);
    }

  } else {

    fprintf(stderr,"Error in gpuBCSlipWallU(): not implemented for ndims=%d.\n", ndims);
    return 1;

  }


  return 0;
}

#endif
