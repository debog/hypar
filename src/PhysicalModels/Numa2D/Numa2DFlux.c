#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/numa2d.h>
#include <hypar.h>

int Numa2DFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar  *solver = (HyPar*)   s;
  Numa2D *param  = (Numa2D*) solver->physics;
  int     i;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
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

    _GetCoordinate_(_YDIR_,index[_YDIR_]-ghosts,dim,ghosts,solver->x,ycoord);
    param->StandardAtmosphere(param,ycoord,&EP,&P0,&rho0,&T0);

    _Numa2DGetFlowVars_     ((u+_MODEL_NVARS_*p),drho,uvel,vvel,dT,rho0);
    _Numa2DComputePressure_ (param,T0,dT,P0,dP);
    _Numa2DSetFlux_         ((f+_MODEL_NVARS_*p),dir,drho,uvel,vvel,dT,dP,rho0,T0);

    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}

int Numa2DStiffFlux(double *f,double *u,int dir,void *s,double t)
{
  HyPar  *solver = (HyPar*)   s;
  Numa2D *param  = (Numa2D*) solver->physics;
  int     i;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
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

    _GetCoordinate_(_YDIR_,index[_YDIR_]-ghosts,dim,ghosts,solver->x,ycoord);
    param->StandardAtmosphere(param,ycoord,&EP,&P0,&rho0,&T0);

    _Numa2DGetFlowVars_               ((u+_MODEL_NVARS_*p),drho,uvel,vvel,dT,rho0);
    _Numa2DComputeLinearizedPressure_ (param,T0,dT,P0,dP);
    _Numa2DSetLinearFlux_             ((f+_MODEL_NVARS_*p),dir,drho,uvel,vvel,dT,dP,rho0,T0);

    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
