/*! @file BCSWSlipWall.c
    @author Debojyoti Ghosh
    @brief Slip-wall boundary conditions for shallow water equations
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

#include <physicalmodels/shallowwater1d.h>
#include <physicalmodels/shallowwater2d.h>

/*! Applies the slip-wall boundary condition: This is specific to the one and two dimenstional
    shallow water equations (#ShallowWater1D, #ShallowWater2D).
    It is used for simulating inviscid walls or symmetric boundaries. The height,
    and tangential velocity at the ghost points are extrapolated from the interior, while the
    normal velocity at the ghost points is set such that the interpolated value at the boundary
    face is equal to the specified wall velocity.
*/
int BCSWSlipWallU(
                 void    *a_b,     /*!< Boundary object of type #DomainBoundary */
                 void    *a_m,     /*!< MPI object of type #MPIVariables */
                 int     a_ndims,  /*!< Number of spatial dimensions */
                 int     a_nvars,  /*!< Number of variables/DoFs per grid point */
                 int     *a_size,  /*!< Integer array with the number of grid points in each spatial dimension */
                 int     a_ghosts, /*!< Number of ghost points */
                 double  *a_phi,   /*!< The solution array on which to apply the boundary condition */
                 double  a_waqt    /*!< Current solution time */
               )
{
  DomainBoundary *boundary = (DomainBoundary*) a_b;

  int dim   = boundary->m_dim;
  int face  = boundary->m_face;

  if (a_ndims == 1) {

    if (boundary->m_on_this_proc) {
      int bounds[a_ndims], indexb[a_ndims], indexi[a_ndims];
      _ArraySubtract1D_(bounds,boundary->m_ie,boundary->m_is,a_ndims);
      _ArraySetValue_(indexb,a_ndims,0);
      int done = 0;
      while (!done) {
        int p1, p2;
        _ArrayCopy1D_(indexb,indexi,a_ndims);
        _ArrayAdd1D_(indexi,indexi,boundary->m_is,a_ndims);
        if      (face ==  1) indexi[dim] = a_ghosts-1-indexb[dim];
        else if (face == -1) indexi[dim] = a_size[dim]-indexb[dim]-1;
        else return(1);
        _ArrayIndex1DWO_(a_ndims,a_size,indexb,boundary->m_is,a_ghosts,p1);
        _ArrayIndex1D_(a_ndims,a_size,indexi,a_ghosts,p2);

        /* flow variables in the interior */
        double h, uvel;
        double h_gpt, uvel_gpt;
        _ShallowWater1DGetFlowVar_((a_phi+a_nvars*p2),h,uvel);
        /* set the ghost point values */
        h_gpt = h;
        uvel_gpt = 2.0*boundary->m_FlowVelocity[_XDIR_] - uvel;

        a_phi[a_nvars*p1+0] = h_gpt;
        a_phi[a_nvars*p1+1] = h_gpt * uvel_gpt;

        _ArrayIncrementIndex_(a_ndims,bounds,indexb,done);
      }
    }

  } else if (a_ndims == 2) {

    if (boundary->m_on_this_proc) {
      int bounds[a_ndims], indexb[a_ndims], indexi[a_ndims];
      _ArraySubtract1D_(bounds,boundary->m_ie,boundary->m_is,a_ndims);
      _ArraySetValue_(indexb,a_ndims,0);
      int done = 0;
      while (!done) {
        int p1, p2;
        _ArrayCopy1D_(indexb,indexi,a_ndims);
        _ArrayAdd1D_(indexi,indexi,boundary->m_is,a_ndims);
        if      (face ==  1) indexi[dim] = a_ghosts-1-indexb[dim];
        else if (face == -1) indexi[dim] = a_size[dim]-indexb[dim]-1;
        else return(1);
        _ArrayIndex1DWO_(a_ndims,a_size,indexb,boundary->m_is,a_ghosts,p1);
        _ArrayIndex1D_(a_ndims,a_size,indexi,a_ghosts,p2);

        /* flow variables in the interior */
        double h, uvel, vvel;
        double h_gpt, uvel_gpt, vvel_gpt;
        _ShallowWater2DGetFlowVar_((a_phi+a_nvars*p2),h,uvel,vvel);
        /* set the ghost point values */
        h_gpt = h;
        if (dim == _XDIR_) {
          uvel_gpt = 2.0*boundary->m_FlowVelocity[_XDIR_] - uvel;
          vvel_gpt = vvel;
        } else if (dim == _YDIR_) {
          uvel_gpt = uvel;
          vvel_gpt = 2.0*boundary->m_FlowVelocity[_YDIR_] - vvel;
        } else {
          uvel_gpt = 0.0;
          vvel_gpt = 0.0;
        }

        a_phi[a_nvars*p1+0] = h_gpt;
        a_phi[a_nvars*p1+1] = h_gpt * uvel_gpt;
        a_phi[a_nvars*p1+2] = h_gpt * vvel_gpt;

        _ArrayIncrementIndex_(a_ndims,bounds,indexb,done);
      }
    }

  }
  return(0);
}
