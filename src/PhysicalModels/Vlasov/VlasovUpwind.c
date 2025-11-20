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

double VlasovAdvectionCoeff(int*, int, void*);

/*! Upwinding scheme for the Vlasov equations */
int VlasovUpwind(  double* a_fI,   /*!< Computed upwind interface flux */
                   double* a_fL,   /*!< Left-biased reconstructed interface flux */
                   double* a_fR,   /*!< Right-biased reconstructed interface flux */
                   double* a_uL,   /*!< Left-biased reconstructed interface solution */
                   double* a_uR,   /*!< Right-biased reconstructed interface solution */
                   double* a_u,    /*!< Cell-centered solution */
                   int     a_dir,  /*!< Spatial dimension */
                   void*   a_s,    /*!< Solver object of type #HyPar */
                   double  a_t   /*!< Current solution time */
                )
{
  HyPar  *solver = (HyPar*)  a_s;
  Vlasov *param  = (Vlasov*) solver->m_physics;
  int    done;

  int ndims   = solver->m_ndims;
  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;

  int index_outer[ndims],
      index_inter[ndims],
      bounds_outer[ndims],
      bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[a_dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[a_dir] += 1;

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {

      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      int indexL[ndims]; _ArrayCopy1D_(index_inter,indexL,ndims); indexL[a_dir]--;
      int indexR[ndims]; _ArrayCopy1D_(index_inter,indexR,ndims);
      int idxL[ndims], idxR[ndims];

      int neig = 2+4*(ndims-1);
      double eig[neig];
      int count = 0;
      eig[count] = VlasovAdvectionCoeff(indexL, a_dir, solver); count++;
      eig[count] = VlasovAdvectionCoeff(indexR, a_dir, solver); count++;
      for (int tdir = 0; tdir < ndims; tdir++ ) {
        if (tdir != a_dir) {
          _ArrayCopy1D_(indexL, idxL, ndims); idxL[tdir]--;
          _ArrayCopy1D_(indexR, idxR, ndims); idxR[tdir]--;
          eig[count] = VlasovAdvectionCoeff(idxL, a_dir, solver); count++;
          eig[count] = VlasovAdvectionCoeff(idxR, a_dir, solver); count++;
          _ArrayCopy1D_(indexL, idxL, ndims); idxL[tdir]++;
          _ArrayCopy1D_(indexR, idxR, ndims); idxR[tdir]++;
          eig[count] = VlasovAdvectionCoeff(idxL, a_dir, solver); count++;
          eig[count] = VlasovAdvectionCoeff(idxR, a_dir, solver); count++;
        }
      }
      if (count != neig) {
        fprintf(stderr, "Error in VlasovUpwind(): count != neig for a_dir=%d\n",a_dir);
        return 1;
      }

      int all_positive = 1, all_negative = 1;
      double maxabs_eig = 0.0;
      for (int n = 0; n<neig; n++) {
        if (eig[n] <= 0) all_positive = 0;
        if (eig[n] >= 0) all_negative = 0;
        if (absolute(eig[n]) > maxabs_eig) {
          maxabs_eig = absolute(eig[n]);
        }
      }

      if (all_positive) {
        a_fI[p] = a_fL[p];
      } else if (all_negative) {
        a_fI[p] = a_fR[p];
      } else {
        a_fI[p] = 0.5 * (a_fL[p] + a_fR[p] - maxabs_eig * (a_uR[p] - a_uL[p]));
      }

    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);

  }

  return 0;
}
