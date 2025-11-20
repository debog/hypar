#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fpdoublewell.h>
#include <hypar.h>

int FPDoubleWellUpwind(double *a_fI,double *a_fL,double *a_fR,double *a_uL,double *a_uR,
                       double *a_u,int a_dir,void *a_s,double a_t)
{
  HyPar        *solver = (HyPar*) a_s;
  int          done,v;

  int ndims   = solver->m_ndims;
  int nvars   = solver->m_nvars;
  int ghosts  = solver->m_ghosts;
  int *dim    = solver->m_dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[a_dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[a_dir] += 1;

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      double x1; _GetCoordinate_(0,index_inter[0]-1,dim,ghosts,solver->m_x,x1);
      double x2; _GetCoordinate_(0,index_inter[0]  ,dim,ghosts,solver->m_x,x2);
      double x = 0.5 * ( x1 + x2 );
      for (v = 0; v < nvars; v++)
        a_fI[nvars*p+v] = (drift(x) > 0 ? a_fL[nvars*p+v] : a_fR[nvars*p+v] );
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
