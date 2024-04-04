/*! @file LinearADRUpwind.c
    @author Debojyoti Ghosh
    @brief Evaluate the upwind flux for the
           linear advection-diffusion-reaction model
*/

#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/linearadr.h>
#include <hypar.h>

/*! Upwinding scheme for linear advection */
int LinearADRUpwind(  double  *fI,  /*!< Computed upwind interface flux */
                      double  *fL,  /*!< Left-biased reconstructed interface flux */
                      double  *fR,  /*!< Right-biased reconstructed interface flux */
                      double  *uL,  /*!< Left-biased reconstructed interface solution */
                      double  *uR,  /*!< Right-biased reconstructed interface solution */
                      double  *u,   /*!< Cell-centered solution */
                      int     dir,  /*!< Spatial dimension */
                      void    *s,   /*!< Solver object of type #HyPar */
                      double  t     /*!< Current solution time */
                   )
{
  HyPar     *solver = (HyPar*)      s;
  LinearADR *param  = (LinearADR*)  solver->physics;
  int       done,v;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int ghosts= solver->ghosts;
  int *dim  = solver->dim_local;

  double *a = param->a;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  if (param->constant_advection == 1) {

    done = 0; _ArraySetValue_(index_outer,ndims,0);
    while (!done) {
      _ArrayCopy1D_(index_outer,index_inter,ndims);
      for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
        int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
        for (v = 0; v < nvars; v++) {
          fI[nvars*p+v] = (a[nvars*dir+v] > 0 ? fL[nvars*p+v] : fR[nvars*p+v] );
        }
      }
      _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
    }

  } else if (param->constant_advection == 0) {

    done = 0; _ArraySetValue_(index_outer,ndims,0);
    while (!done) {
      _ArrayCopy1D_(index_outer,index_inter,ndims);
      for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
        int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
        int indexL[ndims]; _ArrayCopy1D_(index_inter,indexL,ndims); indexL[dir]--;
        int indexR[ndims]; _ArrayCopy1D_(index_inter,indexR,ndims);
        int pL; _ArrayIndex1D_(ndims,dim,indexL,ghosts,pL);
        int pR; _ArrayIndex1D_(ndims,dim,indexR,ghosts,pR);
        for (v = 0; v < nvars; v++) {
          double eigL = a[nvars*ndims*pL+nvars*dir+v],
                 eigR = a[nvars*ndims*pR+nvars*dir+v];
          if ((eigL > 0) && (eigR > 0)) {
             fI[nvars*p+v] = fL[nvars*p+v];
          } else if ((eigL < 0) && (eigR < 0)) {
             fI[nvars*p+v] = fR[nvars*p+v];
          } else {
             double alpha = max(absolute(eigL), absolute(eigR));
             fI[nvars*p+v] = 0.5 * (fL[nvars*p+v] + fR[nvars*p+v] - alpha * (uR[nvars*p+v] - uL[nvars*p+v]));
          }
        }
      }
      _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
    }

  } else {

    done = 0; _ArraySetValue_(index_outer,ndims,0);
    while (!done) {
      _ArrayCopy1D_(index_outer,index_inter,ndims);
      for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
        int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
        for (v = 0; v < nvars; v++) {
          fI[nvars*p+v] = 0.0;
        }
      }
      _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
    }

  }

  return(0);
}

/*! Centered scheme for linear advection */
int LinearADRCenteredFlux(  double  *fI,  /*!< Computed upwind interface flux */
                            double  *fL,  /*!< Left-biased reconstructed interface flux */
                            double  *fR,  /*!< Right-biased reconstructed interface flux */
                            double  *uL,  /*!< Left-biased reconstructed interface solution */
                            double  *uR,  /*!< Right-biased reconstructed interface solution */
                            double  *u,   /*!< Cell-centered solution */
                            int     dir,  /*!< Spatial dimension */
                            void    *s,   /*!< Solver object of type #HyPar */
                            double  t     /*!< Current solution time */
                         )
{
  HyPar     *solver = (HyPar*)      s;
  LinearADR *param  = (LinearADR*)  solver->physics;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int ghosts= solver->ghosts;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  int done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      for (int v = 0; v < nvars; v++) {
        fI[nvars*p+v] = 0.5 * (fL[nvars*p+v] + fR[nvars*p+v]);
      }
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
