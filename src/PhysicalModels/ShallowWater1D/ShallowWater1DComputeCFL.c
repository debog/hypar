/*! @file ShallowWater1DComputeCFL.c
    @author Debojyoti Ghosh
    @brief Contains the function to compute maximum CFL over the domain for the 1D shallow water equations physical model.
*/

#include <stdlib.h>
#include <basic.h>
#include <math.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/shallowwater1d.h>
#include <hypar.h>

/*! Computes the maximum CFL number over the domain. Note that the CFL
    is computed over the local domain on this processor only.
*/
double ShallowWater1DComputeCFL(
                                  void    *a_s, /*!< Solver object of type #HyPar */
                                  void    *a_m, /*!< MPI object of type #MPIVariables */
                                  double  a_dt, /*!< Time step size for which to compute the CFL */
                                  double  a_t   /*!< Time */
                               )
{
  HyPar           *solver = (HyPar*)   a_s;
  ShallowWater1D  *param  = (ShallowWater1D*) solver->m_physics;

  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;
  int ndims   = solver->m_ndims;
  int index[ndims];
  double *u   = solver->m_u;

  double max_cfl = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double h, v, c, dxinv, local_cfl;
    _ShallowWater1DGetFlowVar_((u+_MODEL_NVARS_*p),h,v);
    _GetCoordinate_(0,index[0],dim,ghosts,solver->m_dxinv,dxinv); /* 1/dx */
    c = sqrt(param->m_g*h); /* speed of gravity waves */
    local_cfl = (absolute(v)+c)*a_dt*dxinv; /* local cfl for this grid point */
    if (local_cfl > max_cfl) max_cfl = local_cfl;
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
