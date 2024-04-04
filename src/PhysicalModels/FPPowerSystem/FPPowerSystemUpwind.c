#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem.h>
#include <hypar.h>

double FPPowerSystemDriftFunction(int,void*,double,double,double);

int FPPowerSystemUpwind(double *fI,double *fL,double *fR,double *uL,double *uR,
                        double *u,int dir,void *s,double t)
{
  HyPar         *solver = (HyPar*) s;
  FPPowerSystem *params = (FPPowerSystem*)solver->physics;
  int           done,v;

  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int ghosts  = solver->ghosts;
  int *dim    = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0, p);
      double x = 0,y = 0; /* x,y coordinates of the interface */
      if (dir == 0) {
        /* x-interface */
        double x1, x2;
        _GetCoordinate_(0,index_inter[0]-1,dim,ghosts,solver->x,x1);
        _GetCoordinate_(0,index_inter[0]  ,dim,ghosts,solver->x,x2);
        x = 0.5 * ( x1 + x2 );
        _GetCoordinate_(1,index_inter[1],dim,ghosts,solver->x,y);
      } else if (dir == 1) {
        /* y-interface */
        _GetCoordinate_(0,index_inter[0],dim,ghosts,solver->x,x);
        double y1, y2;
        _GetCoordinate_(1,index_inter[1]-1,dim,ghosts,solver->x,y1);
        _GetCoordinate_(1,index_inter[1]  ,dim,ghosts,solver->x,y2);
        y = 0.5 * ( y1 + y2 );
      }
      double drift = FPPowerSystemDriftFunction(dir,params,x,y,t);
      for (v = 0; v < nvars; v++)
        fI[nvars*p+v] = drift * (drift > 0 ? fL[nvars*p+v] : fR[nvars*p+v] );
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
