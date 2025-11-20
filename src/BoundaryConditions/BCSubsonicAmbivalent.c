/*! @file BCSubsonicAmbivalent.c
    @author Debojyoti Ghosh
    @brief Subsonic "ambivalent" boundary conditions (specific to Euler/Navier-Stokes systems)
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

#include <physicalmodels/euler2d.h>
#include <physicalmodels/navierstokes3d.h>

/*! Applies the subsonic "ambivalent" boundary condition:
    The velocity at the boundary face is extrapolated from the interior and
    its dot product with the boundary normal (pointing into the domain) is
    computed.
    + If it is positive, subsonic inflow boundary conditions are applied: the
      density and velocity at the physical boundary ghost points are specified,
      while the pressure is extrapolated from the interior of the domain.
    + If it is negative, subsonic outflow boundary conditions are applied: the
      pressure at the physical boundary ghost points is specified, while the
      density and velocity are extrapolated from the interior.

    This boundary condition is specific to two and three dimension Euler and
    Navier-Stokes systems (#Euler2D, #NavierStokes2D, #NavierStokes3D).
*/
int BCSubsonicAmbivalentU(
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

  if (a_ndims == 2) {

    /* create a fake physics object */
    Euler2D physics;
    double gamma;
    gamma = physics.m_gamma = boundary->m_gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);

    /* boundary normal (pointing into the domain) */
    double nx, ny;
    if (dim == 0) {
      nx = 1.0;
      ny = 0.0;
    } else if (dim == 1) {
      nx = 0.0;
      ny = 1.0;
    }
    nx *= (double) face;
    ny *= (double) face;

    if (boundary->m_on_this_proc) {
      int bounds[a_ndims], indexb[a_ndims], indexi[a_ndims], indexj[a_ndims];
      _ArraySubtract1D_(bounds,boundary->m_ie,boundary->m_is,a_ndims);
      _ArraySetValue_(indexb,a_ndims,0);
      int done = 0;
      while (!done) {
        int p1, p2;
        double rho, uvel, vvel, energy, pressure;

        /* compute boundary face velocity  - 2nd order */
        _ArrayCopy1D_(indexb,indexi,a_ndims);
        _ArrayAdd1D_(indexi,indexi,boundary->m_is,a_ndims);
        _ArrayCopy1D_(indexi,indexj,a_ndims);
        if (face ==  1) {
          indexi[dim] = 0;
          indexj[dim] = indexi[dim] + 1;
        } else if (face == -1) {
          indexi[dim] = a_size[dim]-1;
          indexj[dim] = indexi[dim] - 1;
        }
        _ArrayIndex1D_(a_ndims,a_size,indexi,a_ghosts,p1);
        _ArrayIndex1D_(a_ndims,a_size,indexj,a_ghosts,p2);
        double uvel1, uvel2, uvelb,
               vvel1, vvel2, vvelb;
        _Euler2DGetFlowVar_((a_phi+a_nvars*p1),rho,uvel1,vvel1,energy,pressure,(&physics));
        _Euler2DGetFlowVar_((a_phi+a_nvars*p2),rho,uvel2,vvel2,energy,pressure,(&physics));
        uvelb = 1.5*uvel1 - 0.5*uvel2;
        vvelb = 1.5*vvel1 - 0.5*vvel2;
        double vel_normal = uvelb*nx + vvelb*ny;

        _ArrayCopy1D_(indexb,indexi,a_ndims);
        _ArrayAdd1D_(indexi,indexi,boundary->m_is,a_ndims);
        if      (face ==  1) indexi[dim] = a_ghosts-1-indexb[dim];
        else if (face == -1) indexi[dim] = a_size[dim]-indexb[dim]-1;
        else return(1);
        _ArrayIndex1DWO_(a_ndims,a_size,indexb,boundary->m_is,a_ghosts,p1);
        _ArrayIndex1D_(a_ndims,a_size,indexi,a_ghosts,p2);

        /* flow variables in the interior */
        _Euler2DGetFlowVar_((a_phi+a_nvars*p2),rho,uvel,vvel,energy,pressure,(&physics));

        /* set the ghost point values */
        double rho_gpt, uvel_gpt, vvel_gpt, energy_gpt, pressure_gpt;
        if (vel_normal > 0) {
          /* inflow */
          rho_gpt = boundary->m_FlowDensity;
          pressure_gpt = pressure;
          uvel_gpt = boundary->m_FlowVelocity[0];
          vvel_gpt = boundary->m_FlowVelocity[1];
        } else {
          /* outflow */
          rho_gpt = rho;
          pressure_gpt = boundary->m_FlowPressure;
          uvel_gpt = uvel;
          vvel_gpt = vvel;
        }
        energy_gpt = inv_gamma_m1*pressure_gpt
                    + 0.5 * rho_gpt * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt);

        a_phi[a_nvars*p1+0] = rho_gpt;
        a_phi[a_nvars*p1+1] = rho_gpt * uvel_gpt;
        a_phi[a_nvars*p1+2] = rho_gpt * vvel_gpt;
        a_phi[a_nvars*p1+3] = energy_gpt;

        _ArrayIncrementIndex_(a_ndims,bounds,indexb,done);
      }
    }

  } else if (a_ndims == 3) {

    /* create a fake physics object */
    double gamma;
    gamma = boundary->m_gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);

    /* boundary normal (pointing into the domain) */
    double nx, ny, nz;
    if (dim == 0) {
      nx = 1.0;
      ny = 0.0;
      nz = 0.0;
    } else if (dim == 1) {
      nx = 0.0;
      ny = 1.0;
      nz = 0.0;
    } else if (dim == 2) {
      nx = 0.0;
      ny = 0.0;
      nz = 1.0;
    }
    nx *= (double) face;
    ny *= (double) face;
    nz *= (double) face;

    if (boundary->m_on_this_proc) {
      int bounds[a_ndims], indexb[a_ndims], indexi[a_ndims], indexj[a_ndims];
      _ArraySubtract1D_(bounds,boundary->m_ie,boundary->m_is,a_ndims);
      _ArraySetValue_(indexb,a_ndims,0);
      int done = 0;
      while (!done) {
        int p1, p2;
        double rho, uvel, vvel, wvel, energy, pressure;

        /* compute boundary face velocity  - 2nd order */
        _ArrayCopy1D_(indexb,indexi,a_ndims);
        _ArrayAdd1D_(indexi,indexi,boundary->m_is,a_ndims);
        _ArrayCopy1D_(indexi,indexj,a_ndims);
        if (face ==  1) {
          indexi[dim] = 0;
          indexj[dim] = indexi[dim] + 1;
        } else if (face == -1) {
          indexi[dim] = a_size[dim]-1;
          indexj[dim] = indexi[dim] - 1;
        }
        _ArrayIndex1D_(a_ndims,a_size,indexi,a_ghosts,p1);
        _ArrayIndex1D_(a_ndims,a_size,indexj,a_ghosts,p2);
        double uvel1, uvel2, uvelb,
               vvel1, vvel2, vvelb,
               wvel1, wvel2, wvelb;
        _NavierStokes3DGetFlowVar_((a_phi+a_nvars*p1),_NavierStokes3D_stride_,rho,uvel1,vvel1,wvel1,energy,pressure,gamma);
        _NavierStokes3DGetFlowVar_((a_phi+a_nvars*p2),_NavierStokes3D_stride_,rho,uvel2,vvel2,wvel2,energy,pressure,gamma);
        uvelb = 1.5*uvel1 - 0.5*uvel2;
        vvelb = 1.5*vvel1 - 0.5*vvel2;
        wvelb = 1.5*wvel1 - 0.5*wvel2;
        double vel_normal = uvelb*nx + vvelb*ny + wvelb*nz;

        _ArrayCopy1D_(indexb,indexi,a_ndims);
        _ArrayAdd1D_(indexi,indexi,boundary->m_is,a_ndims);
        if      (face ==  1) indexi[dim] = a_ghosts-1-indexb[dim];
        else if (face == -1) indexi[dim] = a_size[dim]-indexb[dim]-1;
        else return(1);
        _ArrayIndex1DWO_(a_ndims,a_size,indexb,boundary->m_is,a_ghosts,p1);
        _ArrayIndex1D_(a_ndims,a_size,indexi,a_ghosts,p2);

        /* flow variables in the interior */
        _NavierStokes3DGetFlowVar_((a_phi+a_nvars*p2),_NavierStokes3D_stride_,rho,uvel,vvel,wvel,energy,pressure,gamma);

        /* set the ghost point values */
        double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt;
        if (vel_normal > 0) {
          /* inflow */
          rho_gpt = boundary->m_FlowDensity;
          pressure_gpt = pressure;
          uvel_gpt = boundary->m_FlowVelocity[0];
          vvel_gpt = boundary->m_FlowVelocity[1];
          wvel_gpt = boundary->m_FlowVelocity[2];
        } else {
          /* outflow */
          rho_gpt = rho;
          pressure_gpt = boundary->m_FlowPressure;
          uvel_gpt = uvel;
          vvel_gpt = vvel;
          wvel_gpt = wvel;
        }
        energy_gpt = inv_gamma_m1*pressure_gpt
                    + 0.5 * rho_gpt
                    * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt + wvel_gpt*wvel_gpt);

        a_phi[a_nvars*p1+0] = rho_gpt;
        a_phi[a_nvars*p1+1] = rho_gpt * uvel_gpt;
        a_phi[a_nvars*p1+2] = rho_gpt * vvel_gpt;
        a_phi[a_nvars*p1+3] = rho_gpt * wvel_gpt;
        a_phi[a_nvars*p1+4] = energy_gpt;

        _ArrayIncrementIndex_(a_ndims,bounds,indexb,done);
      }
    }

  }
  return(0);
}
