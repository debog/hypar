/*! @file BCThermalSlipWall.c
    @author Debojyoti Ghosh
    @brief Thermal slip-wall boundary conditions
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>
#include <mpivars.h>

#include <physicalmodels/euler1d.h>
#include <physicalmodels/euler2d.h>
#include <physicalmodels/navierstokes3d.h>

/*! Applies the thermal slip-wall boundary condition: This is specific to the 3D
    Navier-Stokes system (#NavierStokes3D).
    It is used for simulating inviscid walls or symmetric boundaries, where the temperature is
    specified. The  density and the tangential velocity at the ghost points are extrapolated
    for the interior. The normal velocity at the ghost points is set such that the interpolated
    velocity at the boundary face is the specified wall velocity. The pressure at the ghost
    points is set by multiplying the extrapolated density by the specified temperature.

    \b Note: It is assumed that the temperature already contains the gas constant factor, i.e.,
    \f$ T = P/\rho\f$.
*/
int BCThermalSlipWallU(
                         void    *b,     /*!< Boundary object of type #DomainBoundary */
                         void    *m,     /*!< MPI object of type #MPIVariables */
                         int     ndims,  /*!< Number of spatial dimensions */
                         int     nvars,  /*!< Number of variables/DoFs per grid point */
                         int     *size,  /*!< Integer array with the number of grid points in each spatial dimension */
                         int     ghosts, /*!< Number of ghost points */
                         double  *phi,   /*!< The solution array on which to apply the boundary condition */
                         double  waqt    /*!< Current solution time */
                      )
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;
  int face  = boundary->face;

  if (ndims == 3) {

    /* create a fake physics object */
    NavierStokes3D physics; 
    double gamma; 
    gamma = physics.gamma = boundary->gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);

    if (boundary->on_this_proc) {

      int     *temperature_field_size = boundary->UnsteadyTemperatureSize;
      int     n_time_levels           = temperature_field_size[dim];
      double  *time_levels            = boundary->UnsteadyTimeLevels;
      double  *temperature_data       = boundary->UnsteadyTemperatureData;
    
      int it = n_time_levels - 1;
      while ((time_levels[it] > waqt) && (it > 0))  it--;

      int bounds[ndims], indexb[ndims], indexi[ndims];
      _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
      _ArraySetValue_(indexb,ndims,0);

      int done = 0;
      while (!done) {

        int p1, p2;
        _ArrayCopy1D_(indexb,indexi,ndims);
        _ArrayAdd1D_(indexi,indexi,boundary->is,ndims);
        if      (face ==  1) indexi[dim] = ghosts-1-indexb[dim];
        else if (face == -1) indexi[dim] = size[dim]-indexb[dim]-1;
        else return(1);
        _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p1);
        _ArrayIndex1D_(ndims,size,indexi,ghosts,p2);

        /* get the specified temperature */
        int index1[ndims]; _ArrayCopy1D_(indexb,index1,ndims);
        index1[dim] = it;
        int q; _ArrayIndex1D_(ndims,temperature_field_size,index1,0,q);
        double temperature_b = temperature_data[q];
        
        /* flow variables in the interior */
        double rho, uvel, vvel, wvel, energy, pressure;
        _NavierStokes3DGetFlowVar_((phi+nvars*p2),rho,uvel,vvel,wvel,energy,pressure,(&physics));
        /* set the ghost point values */
        double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt;
        rho_gpt = rho;
        if (dim == _XDIR_) {
          uvel_gpt = 2.0*boundary->FlowVelocity[_XDIR_] - uvel;
          vvel_gpt = vvel;
          wvel_gpt = wvel;
        } else if (dim == _YDIR_) {
          uvel_gpt = uvel;
          vvel_gpt = 2.0*boundary->FlowVelocity[_YDIR_] - vvel;
          wvel_gpt = wvel;
        } else if (dim == _ZDIR_) {
          uvel_gpt = uvel;
          vvel_gpt = vvel;
          wvel_gpt = 2.0*boundary->FlowVelocity[_ZDIR_] - wvel;
        } else {
          uvel_gpt = 0.0;
          vvel_gpt = 0.0;
          wvel_gpt = 0.0;
        }
        pressure_gpt = rho_gpt * temperature_b;
        energy_gpt = inv_gamma_m1*pressure_gpt 
                    + 0.5 * rho_gpt 
                    * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt + wvel_gpt*wvel_gpt);

        phi[nvars*p1+0] = rho_gpt;
        phi[nvars*p1+1] = rho_gpt * uvel_gpt;
        phi[nvars*p1+2] = rho_gpt * vvel_gpt;
        phi[nvars*p1+3] = rho_gpt * wvel_gpt;
        phi[nvars*p1+4] = energy_gpt;

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  }
  return(0);
}

/*! Applies the slip-wall boundary condition to the delta-solution: This is specific to the two 
    and three dimensional Euler and Navier-Stokes systems (#Euler2D, #NavierStokes2D, #NavierStokes3D).
    It is used for simulating inviscid walls or symmetric boundaries. The pressure, density,
    and tangential velocity at the ghost points are extrapolated from the interior, while the
    normal velocity at the ghost points is set such that the interpolated value at the boundary 
    face is equal to the specified wall velocity.
    \n\n
    The above treatment is applied on the delta-solution added to the reference solution.
*/
int BCThermalSlipWallDU(
                         void    *b,       /*!< Boundary object of type #DomainBoundary */
                         void    *m,       /*!< MPI object of type #MPIVariables */
                         int     ndims,    /*!< Number of spatial dimensions */
                         int     nvars,    /*!< Number of variables/DoFs per grid point */
                         int     *size,    /*!< Integer array with the number of grid points in each spatial dimension */
                         int     ghosts,   /*!< Number of ghost points */
                         double  *phi,     /*!< The solution array on which to apply the boundary condition -
                                                Note that this is a delta-solution \f$\Delta {\bf U}\f$.*/
                         double  *phi_ref, /*!< Reference solution */
                         double  waqt      /*!< Current solution time */
                       )
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  int dim   = boundary->dim;
  int face  = boundary->face;
  int v;

  if (ndims == 3) {

    /* create a fake physics object */
    NavierStokes3D physics; 
    double gamma; 
    gamma = physics.gamma = boundary->gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);

    if (boundary->on_this_proc) {
      int bounds[ndims], indexb[ndims], indexi[ndims];
      _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
      _ArraySetValue_(indexb,ndims,0);
      int done = 0;
      while (!done) {
        int p1, p2;
        _ArrayCopy1D_(indexb,indexi,ndims);
        _ArrayAdd1D_(indexi,indexi,boundary->is,ndims);
        if      (face ==  1) indexi[dim] = ghosts-1-indexb[dim];
        else if (face == -1) indexi[dim] = size[dim]-indexb[dim]-1;
        else return(1);
        _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p1);
        _ArrayIndex1D_(ndims,size,indexi,ghosts,p2);
        
        /* flow in the interior is phi + phi_ref (since phi is DU) */
        double phi_total[nvars]; 
        for (v=0; v<nvars; v++) phi_total[v] = phi[nvars*p2+v]+phi_ref[nvars*p2+v];
        
        /* flow variables in the interior */
        double rho, uvel, vvel, wvel, energy, pressure;
        double rho0, uvel0, vvel0, wvel0, energy0, pressure0;
        _NavierStokes3DGetFlowVar_(phi_total,rho,uvel,vvel,wvel,energy,pressure,(&physics));
        _NavierStokes3DGetFlowVar_((phi_ref+nvars*p2),rho0,uvel0,vvel0,wvel0,energy0,pressure0,(&physics));
        /* set the ghost point values */
        double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt;
        double rho0_gpt, uvel0_gpt, vvel0_gpt, wvel0_gpt, energy0_gpt, pressure0_gpt;
        /* setting the ghost point values for the total flow variables */
        rho_gpt = rho;
        pressure_gpt = pressure;
        if (dim == _XDIR_) {
          uvel_gpt = 2.0*boundary->FlowVelocity[_XDIR_] - uvel;
          vvel_gpt = vvel;
          wvel_gpt = wvel;
        } else if (dim == _YDIR_) {
          uvel_gpt = uvel;
          vvel_gpt = 2.0*boundary->FlowVelocity[_YDIR_] - vvel;
          wvel_gpt = wvel;
        } else if (dim == _ZDIR_) {
          uvel_gpt = uvel;
          vvel_gpt = vvel;
          wvel_gpt = 2.0*boundary->FlowVelocity[_ZDIR_] - wvel;
        } else {
          uvel_gpt = 0.0;
          vvel_gpt = 0.0;
          wvel_gpt = 0.0;
        }
        energy_gpt = inv_gamma_m1*pressure_gpt 
                    + 0.5 * rho_gpt 
                    * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt + wvel_gpt*wvel_gpt);
        /* setting the ghost point values for the reference flow variables */
        rho0_gpt = rho0;
        pressure0_gpt = pressure0;
        if (dim == _XDIR_) {
          uvel0_gpt = 2.0*boundary->FlowVelocity[_XDIR_] - uvel0;
          vvel0_gpt = vvel0;
          wvel0_gpt = wvel0;
        } else if (dim == _YDIR_) {
          uvel0_gpt = uvel0;
          vvel0_gpt = 2.0*boundary->FlowVelocity[_YDIR_] - vvel0;
          wvel0_gpt = wvel0;
        } else if (dim == _ZDIR_) {
          uvel0_gpt = uvel0;
          vvel0_gpt = vvel0;
          wvel0_gpt = 2.0*boundary->FlowVelocity[_ZDIR_] - wvel0;
        } else {
          uvel0_gpt = 0.0;
          vvel0_gpt = 0.0;
          wvel0_gpt = 0.0;
        }
        energy0_gpt = inv_gamma_m1*pressure0_gpt 
                    + 0.5 * rho0_gpt 
                    * (uvel0_gpt*uvel0_gpt + vvel0_gpt*vvel0_gpt + wvel0_gpt*wvel0_gpt);

        phi[nvars*p1+0] = rho_gpt            - rho0_gpt;
        phi[nvars*p1+1] = rho_gpt * uvel_gpt - rho0_gpt * uvel0_gpt;
        phi[nvars*p1+2] = rho_gpt * vvel_gpt - rho0_gpt * vvel0_gpt;
        phi[nvars*p1+3] = rho_gpt * wvel_gpt - rho0_gpt * wvel0_gpt;
        phi[nvars*p1+4] = energy_gpt         - energy0_gpt;

        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }

  }
  return(0);
}
