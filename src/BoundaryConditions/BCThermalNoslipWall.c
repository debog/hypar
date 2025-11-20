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

  if (a_ndims == 3) {

    /* create a fake physics object */
    double gamma;
    gamma = boundary->m_gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);

    if (boundary->m_on_this_proc) {

      int     *temperature_field_size = boundary->m_UnsteadyTemperatureSize;
      int     n_time_levels           = temperature_field_size[dim];
      double  *time_levels            = boundary->m_UnsteadyTimeLevels;
      double  *temperature_data       = boundary->m_UnsteadyTemperatureData;

      int it = n_time_levels - 1;
      while ((time_levels[it] > a_waqt) && (it > 0))  it--;

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

        /* get the specified temperature */
        int index1[a_ndims]; _ArrayCopy1D_(indexb,index1,a_ndims);
        index1[dim] = it;
        int q; _ArrayIndex1D_(a_ndims,temperature_field_size,index1,0,q);
        double temperature_b = temperature_data[q];

        /* flow variables in the interior */
        double rho, uvel, vvel, wvel, energy, pressure;
        _NavierStokes3DGetFlowVar_((a_phi+a_nvars*p2),_NavierStokes3D_stride_,rho,uvel,vvel,wvel,energy,pressure,gamma);
        /* set the ghost point values */
        double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt;
        rho_gpt = rho;
        uvel_gpt = 2.0*boundary->m_FlowVelocity[0] - uvel;
        vvel_gpt = 2.0*boundary->m_FlowVelocity[1] - vvel;
        wvel_gpt = 2.0*boundary->m_FlowVelocity[2] - wvel;
        pressure_gpt = rho_gpt * temperature_b;
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
