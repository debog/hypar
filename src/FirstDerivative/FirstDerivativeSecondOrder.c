#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <secondderivative.h>
#include <mpivars.h>
#include <hypar.h>

/* 
  Second order central differencing
*/

int FirstDerivativeSecondOrderCentral(double *Df,double *f,int dir,void *s,void *m)
{
  HyPar         *solver = (HyPar*) s;
  int           i, v;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;


  if ((!Df) || (!f)) {
    fprintf(stderr, "Error in FirstDerivativeSecondOrder(): input arrays not allocated.\n");
    return(1);
  }

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], index_outer[ndims], bounds_outer[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;

  int done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,indexC,ndims);
    /* left boundary */
    for (i = -ghosts; i < -ghosts+1; i++) {
      int qC, qR, qRR;
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qR );
      indexC[dir] = i+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qRR);
      for (v=0; v<nvars; v++)  Df[qC*nvars+v] = -3*f[qC*nvars+v]+4*f[qR*nvars+v]-f[qRR*nvars+v];
    }
    /* interior */
    for (i = -ghosts+1; i < dim[dir]+ghosts-1; i++) {
      int qC, qL, qR;
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qL);
      indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qR);
      for (v=0; v<nvars; v++)  Df[qC*nvars+v] = f[qR*nvars+v]-f[qL*nvars+v];
    }
    /* right boundary */
    for (i = dim[dir]+ghosts-1; i < dim[dir]+ghosts; i++) {
      int qLL, qL, qC;
      indexC[dir] = i-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qLL);
      indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qL );
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      for (v=0; v<nvars; v++)  Df[qC*nvars+v] = 3*f[qC*nvars+v]-4*f[qL*nvars+v]+f[qLL*nvars+v];
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }
  
  return(0);
}