#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/burgers.h>
#include <hypar.h>

int BurgersUpwind(double *fI,double *fL,double *fR,double *uL,double *uR,
                  double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)   s;
  Burgers   *param  = (Burgers*) solver->physics;
  int        done,v;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      for (v = 0; v < nvars; v++) {
        fI[nvars*p+v] = (u[nvars*p+v] > 0 ? fL[nvars*p+v] : fR[nvars*p+v]);
#if 0
        double uj = u[nvars*p+v];
        if (uj > 0) {
           fI[nvars*p+v] = fL[nvars*p+v];
        }
        else if (uj < 0) {
           fI[nvars*p+v] = fR[nvars*p+v];
        }
        else { // uj == 0
           double ujp = u[nvars*p+v+1];
           double u_max_abs = max(abs(uj), abs(ujp));
           fI[nvars*p+v] = 0.5 * (fL[nvars*p+v] + fR[nvars*p+v]) -
              u_max_abs * (uR[nvars*p+v] - uL[nvars*p+v]);
        }
#endif
      }
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
