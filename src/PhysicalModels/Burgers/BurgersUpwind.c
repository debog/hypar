/*! @file BurgersUpwind.c
    @author John Loffeld
    @brief Contains functions to compute the upwind flux at grid interfaces for the Burgers equations.
*/

#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/burgers.h>
#include <hypar.h>

/*! Upwinding scheme for the Burgers equations */
int BurgersUpwind(  double* a_fI,   /*!< Computed upwind interface flux */
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
  HyPar   *solver = (HyPar*)   a_s;
  Burgers *param  = (Burgers*) solver->m_physics;
  int     done,v;

  int ndims   = solver->m_ndims;
  int nvars   = solver->m_nvars;
  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[a_dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[a_dir] += 1;

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {

      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      int indexL[ndims]; _ArrayCopy1D_(index_inter,indexL,ndims); indexL[a_dir]--;
      int indexR[ndims]; _ArrayCopy1D_(index_inter,indexR,ndims);
      int pL; _ArrayIndex1D_(ndims,dim,indexL,ghosts,pL);
      int pR; _ArrayIndex1D_(ndims,dim,indexR,ghosts,pR);

      for (v = 0; v < nvars; v++) {
        double eigL = a_u[nvars*pL+v],
               eigR = a_u[nvars*pR+v];

        if ((eigL > 0) && (eigR > 0)) {
           a_fI[nvars*p+v] = a_fL[nvars*p+v];
        } else if ((eigL < 0) && (eigR < 0)) {
           a_fI[nvars*p+v] = a_fR[nvars*p+v];
        } else {
           double alpha = max(absolute(eigL), absolute(eigR));
           a_fI[nvars*p+v] = 0.5 * (a_fL[nvars*p+v] + a_fR[nvars*p+v] - alpha * (a_uR[nvars*p+v] - a_uL[nvars*p+v]));
        }
      }

    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);

  }

  return(0);
}
