/*! @file BCThermalNoslipWall.c
    @author Debojyoti Ghosh
    @brief Thermal no-slip-wall boundary conditions
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>
#include <mpivars.h>

#include <physicalmodels/euler1d.h>
#include <physicalmodels/euler2d.h>
#include <physicalmodels/navierstokes3d.h>

/*! Applies the thermal no-slip-wall boundary condition: This is specific to the 3D
    Navier-Stokes system (#NavierStokes3D).
    It is used for simulating walls boundaries, where the temperature is specified.
    The density at the ghost points is extrapolated from the interior. The
    velocity at the ghost points is set such that the interpolated
    velocity at the boundary face is the specified wall velocity. The pressure at the ghost
    points is set by multiplying the extrapolated density by the specified temperature.

    \b Note: It is assumed that the temperature already contains the gas constant factor, i.e.,
    \f$ T = P/\rho\f$.
*/
int BCThermalNoslipWallU(
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
    double gamma;
    gamma = boundary->gamma;
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
        _NavierStokes3DGetFlowVar_((phi+nvars*p2),_NavierStokes3D_stride_,rho,uvel,vvel,wvel,energy,pressure,gamma);
        /* set the ghost point values */
        double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt;
        rho_gpt = rho;
        uvel_gpt = 2.0*boundary->FlowVelocity[0] - uvel;
        vvel_gpt = 2.0*boundary->FlowVelocity[1] - vvel;
        wvel_gpt = 2.0*boundary->FlowVelocity[2] - wvel;
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
