/*! @file BCSupersonicInflow.c
    @author Debojyoti Ghosh
    @brief Supersonic inflow boundary conditions (specific to Euler/Navier-Stokes systems)
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

#include <physicalmodels/euler2d.h>
#include <physicalmodels/navierstokes3d.h>

/*! Applies the supersonic (steady) inflow boundary condition: All the flow variables
    (density, pressure, velocity) are specified at the physical boundary ghost points,
    since it is supersonic inflow. This boundary condition is specific to two and three
    dimensional Euler/Navier-Stokes systems (#Euler2D, #NavierStokes2D, #NavierStokes3D).
    \n\n
    Note: the Dirichlet boundary condition (#_DIRICHLET_) could be used as well for
    supersonic inflow; however the specified Dirichlet state should be in terms of the
    conserved variables, while the specified supersonic inflow state here is in terms of
    the flow variables.
*/
int BCSupersonicInflowU(
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

  if (a_ndims == 2) {

    double gamma;
    gamma = boundary->m_gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);

    if (boundary->m_on_this_proc) {
      int bounds[a_ndims], indexb[a_ndims];
      _ArraySubtract1D_(bounds,boundary->m_ie,boundary->m_is,a_ndims);
      _ArraySetValue_(indexb,a_ndims,0);
      int done = 0;
      while (!done) {
        int p1; _ArrayIndex1DWO_(a_ndims,a_size,indexb,boundary->m_is,a_ghosts,p1);

        /* set the ghost point values */
        double rho_gpt, uvel_gpt, vvel_gpt, energy_gpt, pressure_gpt;
        rho_gpt      = boundary->m_FlowDensity;
        pressure_gpt = boundary->m_FlowPressure;
        uvel_gpt     = boundary->m_FlowVelocity[0];
        vvel_gpt     = boundary->m_FlowVelocity[1];
        energy_gpt   = inv_gamma_m1*pressure_gpt
                       + 0.5 * rho_gpt * (uvel_gpt*uvel_gpt + vvel_gpt*vvel_gpt);

        a_phi[a_nvars*p1+0] = rho_gpt;
        a_phi[a_nvars*p1+1] = rho_gpt * uvel_gpt;
        a_phi[a_nvars*p1+2] = rho_gpt * vvel_gpt;
        a_phi[a_nvars*p1+3] = energy_gpt;

        _ArrayIncrementIndex_(a_ndims,bounds,indexb,done);
      }
    }

  } else if (a_ndims == 3) {

    double gamma;
    gamma = boundary->m_gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);

    if (boundary->m_on_this_proc) {
      int bounds[a_ndims], indexb[a_ndims];
      _ArraySubtract1D_(bounds,boundary->m_ie,boundary->m_is,a_ndims);
      _ArraySetValue_(indexb,a_ndims,0);
      int done = 0;
      while (!done) {
        int p1; _ArrayIndex1DWO_(a_ndims,a_size,indexb,boundary->m_is,a_ghosts,p1);

        /* set the ghost point values */
        double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt;
        rho_gpt      = boundary->m_FlowDensity;
        pressure_gpt = boundary->m_FlowPressure;
        uvel_gpt     = boundary->m_FlowVelocity[0];
        vvel_gpt     = boundary->m_FlowVelocity[1];
        wvel_gpt     = boundary->m_FlowVelocity[2];
        energy_gpt   = inv_gamma_m1*pressure_gpt
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
