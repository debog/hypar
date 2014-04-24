#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <hypar.h>

int FPPowerSystem3BusDissipationFunction(int,int,void*,double,double*);

int FPPowerSystem3BusDiffusion(double *f,double *u,int dir1,int dir2,void *s,double t)
{
  HyPar             *solver = (HyPar*)              s;
  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*)  solver->physics;
  int               i, v;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double dissipation[ndims*ndims]; FPPowerSystem3BusDissipationFunction(dir1,dir2,params,t,dissipation);
    for (v = 0; v < nvars; v++) f[nvars*p+v] = dissipation[dir1*ndims+dir2] * u[nvars*p+v];
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}