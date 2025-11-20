#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/numa3d.h>
#include <hypar.h>

double Numa3DComputeCFL(void *a_s,void *a_m,double a_dt,double a_t)
{
  HyPar  *solver = (HyPar*)  a_s;
  Numa3D *param  = (Numa3D*) solver->m_physics;

  int     *dim    = solver->m_dim_local;
  int ghosts = solver->m_ghosts;
  int     ndims   = solver->m_ndims;
  double  *u      = solver->m_u;
  int     index[ndims];

  double max_cfl = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double drho,uvel,vvel,wvel,dT,rho0,T0,P0,EP,c,zcoord;

    _GetCoordinate_(_ZDIR_,index[_ZDIR_],dim,ghosts,solver->m_x,zcoord);
    param->StandardAtmosphere(param,zcoord,&EP,&P0,&rho0,&T0);
    _Numa3DGetFlowVars_         ((u+_MODEL_NVARS_*p),drho,uvel,vvel,wvel,dT,rho0);
    _Numa3DComputeSpeedofSound_ (param->m_gamma,param->m_R,T0,dT,rho0,drho,EP,c);

    double dxinv, dyinv, dzinv;
    _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,solver->m_dxinv,dxinv); /* 1/dx */
    _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->m_dxinv,dyinv); /* 1/dy */
    _GetCoordinate_(_ZDIR_,index[_ZDIR_],dim,ghosts,solver->m_dxinv,dzinv); /* 1/dz */

    double local_cfl[3];
    local_cfl[_XDIR_] = (absolute(uvel)+c)*a_dt*dxinv; /* local cfl for this grid point (x) */
    local_cfl[_YDIR_] = (absolute(vvel)+c)*a_dt*dyinv; /* local cfl for this grid point (y) */
    local_cfl[_ZDIR_] = (absolute(wvel)+c)*a_dt*dzinv; /* local cfl for this grid point (z) */
    if (local_cfl[_XDIR_] > max_cfl) max_cfl = local_cfl[_XDIR_];
    if (local_cfl[_YDIR_] > max_cfl) max_cfl = local_cfl[_YDIR_];
    if (local_cfl[_ZDIR_] > max_cfl) max_cfl = local_cfl[_ZDIR_];

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
