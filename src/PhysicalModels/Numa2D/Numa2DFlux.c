#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/numa2d.h>
#include <hypar.h>

int Numa2DFlux(double *a_f,double *a_u,int a_dir,void *a_s,double a_t)
{
  HyPar  *solver = (HyPar*)   a_s;
  Numa2D *param  = (Numa2D*) solver->m_physics;
  int     i;

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
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double drho,uvel,vvel,dT,dP,rho0,T0,P0,EP,ycoord;

    _GetCoordinate_(_YDIR_,index[_YDIR_]-ghosts,dim,ghosts,solver->m_x,ycoord);
    param->StandardAtmosphere(param,ycoord,&EP,&P0,&rho0,&T0);

    _Numa2DGetFlowVars_     ((a_u+_MODEL_NVARS_*p),drho,uvel,vvel,dT,rho0);
    _Numa2DComputePressure_ (param,T0,dT,P0,dP);
    _Numa2DSetFlux_         ((a_f+_MODEL_NVARS_*p),a_dir,drho,uvel,vvel,dT,dP,rho0,T0);

    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}

int Numa2DStiffFlux(double *a_f,double *a_u,int a_dir,void *a_s,double a_t)
{
  HyPar  *solver = (HyPar*)   a_s;
  Numa2D *param  = (Numa2D*) solver->m_physics;
  int     i;

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
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double drho,uvel,vvel,dT,dP,rho0,T0,P0,EP,ycoord;

    _GetCoordinate_(_YDIR_,index[_YDIR_]-ghosts,dim,ghosts,solver->m_x,ycoord);
    param->StandardAtmosphere(param,ycoord,&EP,&P0,&rho0,&T0);

    _Numa2DGetFlowVars_               ((a_u+_MODEL_NVARS_*p),drho,uvel,vvel,dT,rho0);
    _Numa2DComputeLinearizedPressure_ (param,T0,dT,P0,dP);
    _Numa2DSetLinearFlux_             ((a_f+_MODEL_NVARS_*p),a_dir,drho,uvel,vvel,dT,dP,rho0,T0);

    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
