/*! @file NavierStokes3DFunctions.c
    @author Debojyoti Ghosh
    @brief Miscellaneous functions for the 3D Navier Stokes equations.
*/
#include <stdio.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

/*! Compute the Roe-averaged state for the 3D Navier Stokes equations. This function
    just calls the macro #_NavierStokes3DRoeAverage_ and is not used by any
    functions within the 3D Navier Stokes module. However, it's necessary to define it
    and provide it to the the solver object (#HyPar) so that it can then send it
    to interpolation functions for a characteristic-based reconstruction.
*/
int NavierStokes3DRoeAverage(
                              double  *a_uavg, /*!< The computed Roe-averaged state */
                              double  *a_uL,   /*!< Left state (conserved variables)*/
                              double  *a_uR,   /*!< Right state (conserved variables)*/
                              void    *a_p     /*!< Object of type #NavierStokes3D with physics-related variables */
                            )
{
  NavierStokes3D *param  = (NavierStokes3D*) a_p;
  _NavierStokes3DRoeAverage_(a_uavg,_NavierStokes3D_stride_,a_uL,a_uR,param->m_gamma);
  return(0);
}

/*! Compute the pressure from the conserved solution on a grid */
int NavierStokes3DComputePressure(  double*             a_P, /*!< Array to hold the computed pressure (same layout as a_u) */
                                    const double* const a_u, /*!< Array with the solution vector */
                                    void*               a_s  /*!< Solver object of type #HyPar */
                                  )
{
  HyPar          *solver = (HyPar*)   a_s;
  NavierStokes3D *param  = (NavierStokes3D*) solver->m_physics;
  int            i;

  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;
  int ndims   = solver->m_ndims;
  int index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int idx; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,idx);
    double rho, vx, vy, vz, e, pressure;
    _NavierStokes3DGetFlowVar_((a_u+_MODEL_NVARS_*idx),_NavierStokes3D_stride_,rho,vx,vy,vz,e,pressure,param->m_gamma);
    a_P[idx] = pressure;
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}

/*! Compute the temperature from the conserved solution on a grid */
int NavierStokes3DComputeTemperature( double*             a_T, /*!< Array to hold the computed pressure (same layout as a_u) */
                                      const double* const a_u, /*!< Array with the solution vector */
                                      void*               a_s  /*!< Solver object of type #HyPar */
                                    )
{
  HyPar          *solver = (HyPar*)   a_s;
  NavierStokes3D *param  = (NavierStokes3D*) solver->m_physics;
  int            i;

  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;
  int ndims   = solver->m_ndims;
  int index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int idx; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,idx);
    double rho, vx, vy, vz, e, pressure;
    _NavierStokes3DGetFlowVar_((a_u+_MODEL_NVARS_*idx),_NavierStokes3D_stride_,rho,vx,vy,vz,e,pressure,param->m_gamma);
    a_T[idx] = pressure/rho;
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
