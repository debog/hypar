#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/numa2d.h>
#include <hypar.h>

double Numa2DComputeCFL(void *s,void *m,double dt,double t)
{
  HyPar  *solver = (HyPar*)  s;
  Numa2D *param  = (Numa2D*) solver->m_physics;

  int     *dim    = solver->m_dim_local;
  int     ghosts  = solver->m_ghosts;
  int     ndims   = solver->m_ndims;
  double  *u      = solver->m_u;
  int     index[ndims];

  double max_cfl = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double drho,uvel,vvel,dT,rho0,T0,P0,EP,c,ycoord;

    _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->m_x,ycoord);
    param->StandardAtmosphere   (param,ycoord,&EP,&P0,&rho0,&T0);
    _Numa2DGetFlowVars_         ((u+_MODEL_NVARS_*p),drho,uvel,vvel,dT,rho0);
    _Numa2DComputeSpeedofSound_ (param->m_gamma,param->m_R,T0,dT,rho0,drho,EP,c);

    double dxinv, dyinv;
    _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,solver->m_dxinv,dxinv); /* 1/dx */
    _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,solver->m_dxinv,dyinv); /* 1/dy */

    double local_cfl[2];
    local_cfl[_XDIR_] = (absolute(uvel)+c)*dt*dxinv; /* local cfl for this grid point (x) */
    local_cfl[_YDIR_] = (absolute(vvel)+c)*dt*dyinv; /* local cfl for this grid point (y) */
    if (local_cfl[_XDIR_] > max_cfl) max_cfl = local_cfl[_XDIR_];
    if (local_cfl[_YDIR_] > max_cfl) max_cfl = local_cfl[_YDIR_];

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
