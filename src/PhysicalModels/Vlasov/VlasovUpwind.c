/*! @file VlasovUpwind.c
    @author John Loffeld
    @brief Contains functions to compute the upwind flux at grid interfaces for the Vlasov equations.
*/

#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/vlasov.h>
#include <hypar.h>

/*! Upwinding scheme for the Vlasov equations */
int VlasovUpwind(  double* fI,   /*!< Computed upwind interface flux */
                   double* fL,   /*!< Left-biased reconstructed interface flux */
                   double* fR,   /*!< Right-biased reconstructed interface flux */
                   double* uL,   /*!< Left-biased reconstructed interface solution */
                   double* uR,   /*!< Right-biased reconstructed interface solution */
                   double* u,    /*!< Cell-centered solution */
                   int     dir,  /*!< Spatial dimension */
                   void*   s,    /*!< Solver object of type #HyPar */
                   double  t     /*!< Current solution time */
                )
{
  HyPar  *solver = (HyPar*)  s;
  Vlasov *param  = (Vlasov*) solver->physics;
  int    done;

  int ndims   = solver->ndims;
  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;

#ifdef fftw
  double *field = param->field;
#endif

  bool self_consistent_electric_field = param->self_consistent_electric_field;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {

      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      int indexL[ndims]; _ArrayCopy1D_(index_inter,indexL,ndims); indexL[dir]--;
      int indexR[ndims]; _ArrayCopy1D_(index_inter,indexR,ndims);

      /* Vlasov equation has only one component */
      double eigL, eigR;
      if (dir == 0) {
        /* Velocity coordinate is the HyPar "velocity" */
        _GetCoordinate_(1,indexL[1],dim,ghosts,solver->x,eigL);
        _GetCoordinate_(1,indexR[1],dim,ghosts,solver->x,eigR);
      } else {
        if (self_consistent_electric_field) {
#ifdef fftw
          /* assumes field has been calculated just prior in VlasovAdvection */
          eigL = field[indexL[0]];
          eigR = field[indexR[0]];
#endif
        } else {
          /* Prescribed electric field is the HyPar "velocity" */
          double xL, xR;
          _GetCoordinate_(0,indexL[0],dim,ghosts,solver->x,xL);
          _GetCoordinate_(0,indexR[0],dim,ghosts,solver->x,xR);
          eigL = 0.1 * cos(xL);
          eigR = 0.1 * cos(xR);
        }
      }

      if ((eigL > 0) && (eigR > 0)) {
        fI[p] = fL[p];
      } else if ((eigL < 0) && (eigR < 0)) {
        fI[p] = fR[p];
      } else { 
        double alpha = max(abs(eigL), abs(eigR));
        fI[p] = 0.5 * (fL[p] + fR[p] - alpha * (uR[p] - uL[p]));
      }

    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);

  }

  return(0);
}
