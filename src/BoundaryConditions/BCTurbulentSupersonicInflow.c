/*! @file BCTurbulentSupersonicInflow.c
    @author Debojyoti Ghosh
    @brief Turbulent supersonic inflow boundary condition (specific to the 3D Navier-Stokes system).
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>
#include <mpivars.h>

/*! Applies the turbulent supersonic inflow boundary condition: The inflow consists
    of a mean supersonic inflow on which turbulent flow fluctuations are added. This
    boundary condition is specific to the 3D Navier-Stokes system (#NavierStokes3D).
    \n\n
    Note: Some parts of the code may be hardcoded for use with the shock-turbulence
    interaction problem (for which this boundary condition was written).
*/
int BCTurbulentSupersonicInflowU(
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

  double *inflow_data = boundary->UnsteadyDirichletData;
  int    *inflow_size = boundary->UnsteadyDirichletSize;

  if (ndims == 3) {

    /* create a fake physics object */
    double gamma;
    gamma = boundary->gamma;
    double inv_gamma_m1 = 1.0/(gamma-1.0);

    if (boundary->on_this_proc) {
      /* the following bit is hardcoded for the inflow data
       * representing fluctuations in a domain of length 2pi */
      double  xt = boundary->FlowVelocity[dim] * waqt;
      int     N  = inflow_size[dim];
      double  L  = 2.0 * (4.0*atan(1.0));
      int     it = ((int) ((xt/L) * ((double)N))) % N;

      int bounds[ndims], indexb[ndims];
      _ArraySubtract1D_(bounds,boundary->ie,boundary->is,ndims);
      _ArraySetValue_(indexb,ndims,0);
      int done = 0;
      while (!done) {
        int p1; _ArrayIndex1DWO_(ndims,size,indexb,boundary->is,ghosts,p1);

        /* set the ghost point values - mean flow */
        double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt;
        rho_gpt      = boundary->FlowDensity;
        pressure_gpt = boundary->FlowPressure;
        uvel_gpt     = boundary->FlowVelocity[0];
        vvel_gpt     = boundary->FlowVelocity[1];
        wvel_gpt     = boundary->FlowVelocity[2];

        /* calculate the turbulent fluctuations */
        double duvel , dvvel , dwvel ;
        int index1[ndims]; _ArrayCopy1D_(indexb,index1,ndims);
        index1[dim] = it;
        int q; _ArrayIndex1D_(ndims,inflow_size,index1,0,q);
        duvel = inflow_data[q*nvars+1];
        dvvel = inflow_data[q*nvars+2];
        dwvel = inflow_data[q*nvars+3];

        /* add the turbulent fluctuations to the velocity field */
        uvel_gpt      += duvel;
        vvel_gpt      += dvvel;
        wvel_gpt      += dwvel;

        /* set the ghost point values */
        energy_gpt   = inv_gamma_m1*pressure_gpt
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
