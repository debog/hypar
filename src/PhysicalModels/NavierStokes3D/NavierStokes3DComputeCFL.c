/*! @file NavierStokes3DComputeCFL.c
    @author Debojyoti Ghosh
    @brief Compute the maximum CFL.
*/
#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

/*! Computes the maximum CFL number over the domain. Note that the CFL
    is computed over the local domain on this processor only.
*/
double NavierStokes3DComputeCFL(
                                void    *s, /*!< Solver object of type #HyPar */
                                void    *m, /*!< MPI object of type #MPIVariables */
                                double  dt, /*!< Time step size for which to compute the CFL */
                                double  t   /*!< Time */
                               )
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes3D    *param  = (NavierStokes3D*) solver->m_physics;

  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;
  int ndims   = solver->m_ndims;
  int index[ndims];
  double *u   = solver->m_u;

  double max_cfl = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double rho, vx, vy, vz, e, P, c, dxinv, dyinv, dzinv, local_cfl[3];
    _NavierStokes3DGetFlowVar_((u+_MODEL_NVARS_*p),_NavierStokes3D_stride_,rho,vx,vy,vz,e,P,param->m_gamma);

    c = sqrt(param->m_gamma*P/rho); /* speed of sound */
    _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,solver->m_dxinv,dxinv); /* 1/dx */
    _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->m_dxinv,dyinv); /* 1/dy */
    _GetCoordinate_(_ZDIR_,index[_ZDIR_],dim,ghosts,solver->m_dxinv,dzinv); /* 1/dz */

    local_cfl[_XDIR_] = (absolute(vx)+c)*dt*dxinv; /* local cfl for this grid point (x) */
    local_cfl[_YDIR_] = (absolute(vy)+c)*dt*dyinv; /* local cfl for this grid point (y) */
    local_cfl[_ZDIR_] = (absolute(vz)+c)*dt*dzinv; /* local cfl for this grid point (z) */
    if (local_cfl[_XDIR_] > max_cfl) max_cfl = local_cfl[_XDIR_];
    if (local_cfl[_YDIR_] > max_cfl) max_cfl = local_cfl[_YDIR_];
    if (local_cfl[_ZDIR_] > max_cfl) max_cfl = local_cfl[_ZDIR_];

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
