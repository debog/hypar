#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/numa2d.h>
#include <hypar.h>

int Numa2DSource(double *a_S,double *a_u,void *a_s,void *a_m,double a_t)
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
    double drho,uvel,vvel,dT,rho0,P0,EP,T0,ycoord;

    _GetCoordinate_(_YDIR_,index[_YDIR_]-ghosts,dim,ghosts,solver->m_x,ycoord);
    param->StandardAtmosphere (param,ycoord,&EP,&P0,&rho0,&T0);
    _Numa2DGetFlowVars_       ((a_u+_MODEL_NVARS_*p),drho,uvel,vvel,dT,rho0);
    _Numa2DSetSource_         ((a_S+_MODEL_NVARS_*p),param,drho);

    /* some useless statements to avoid compiler warnings */
    uvel = dT;
    vvel = uvel;
    dT = vvel;

    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
