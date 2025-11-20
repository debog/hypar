/*! @file TemplateModelUpwind.c
    @author [YOUR NAME]
    @brief Upwinding scheme for the Template Model
*/

#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/template_model.h>
#include <hypar.h>

/*!
 * Compute the upwind interface flux for the Template Model.
 *
 * This function computes the numerical flux at grid interfaces by combining
 * the left- and right-biased reconstructed fluxes using an upwinding scheme.
 *
 * Common upwinding schemes:
 * 1. Simple upwind: use flux from upwind direction based on wave speed
 * 2. Local Lax-Friedrichs (Rusanov): (f_L + f_R)/2 - alpha*(u_R - u_L)/2
 * 3. Roe scheme: exact solution of Riemann problem with linearization
 * 4. HLLC/AUSM: approximate Riemann solvers
*/
int TemplateModelUpwind(
  double  *fI,    /*!< Computed upwind interface flux */
  double  *fL,    /*!< Left-biased reconstructed interface flux */
  double  *fR,    /*!< Right-biased reconstructed interface flux */
  double  *uL,    /*!< Left-biased reconstructed interface solution */
  double  *uR,    /*!< Right-biased reconstructed interface solution */
  double  *u,     /*!< Cell-centered solution (needed for some schemes) */
  int     dir,    /*!< Spatial dimension */
  void    *s,     /*!< Solver object of type #HyPar */
  double  t       /*!< Current solution time */
)
{
  HyPar         *solver  = (HyPar*)         s;
  TemplateModel *physics = (TemplateModel*) solver->m_physics;

  int ndims  = solver->m_ndims;
  int nvars  = solver->m_nvars;
  int *dim   = solver->m_dim_local;
  int ghosts = solver->m_ghosts;

  /*
   * Interface array indexing:
   * - Interface arrays (fI, fL, fR, uL, uR) have one extra point in direction dir
   * - They do NOT include ghost points
   * - Outer loop iterates over dimensions perpendicular to dir
   * - Inner loop iterates over interfaces in direction dir
  */

  int index_outer[ndims], index_inter[ndims];
  int bounds_outer[ndims], bounds_inter[ndims];

  _ArrayCopy1D_(dim,bounds_outer,ndims);
  bounds_outer[dir] = 1;
  _ArrayCopy1D_(dim,bounds_inter,ndims);
  bounds_inter[dir] += 1;

  int done = 0;
  _ArraySetValue_(index_outer,ndims,0);

  while (!done) {

    _ArrayCopy1D_(index_outer,index_inter,ndims);

    for (index_inter[dir]=0; index_inter[dir]<bounds_inter[dir]; index_inter[dir]++) {

      /* Interface index */
      int p;
      _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);

      /* Indices of left and right cells adjacent to interface */
      int indexL[ndims], indexR[ndims];
      _ArrayCopy1D_(index_inter,indexL,ndims);
      indexL[dir]--;
      _ArrayCopy1D_(index_inter,indexR,ndims);

      int pL, pR;
      _ArrayIndex1D_(ndims,dim,indexL,ghosts,pL);
      _ArrayIndex1D_(ndims,dim,indexR,ghosts,pR);

      /* Compute upwind flux for each variable */
      int v;
      for (v = 0; v < nvars; v++) {

        /* [REPLACE] Calculate wave speeds (eigenvalues) */
        /* For scalar equation: wave speed is flux derivative df/du */
        /* For systems: compute eigenvalues of flux Jacobian df/du */

        /* Example 1: Scalar advection with constant speed a */
        /*
        double lambda = physics->param1;
        if (lambda > 0) {
          fI[nvars*p+v] = fL[nvars*p+v];  // upwind from left
        } else {
          fI[nvars*p+v] = fR[nvars*p+v];  // upwind from right
        }
        */

        /* Example 2: Scalar with variable wave speed */
        /*
        double lambdaL = uL[nvars*p+v];  // wave speed at left state
        double lambdaR = uR[nvars*p+v];  // wave speed at right state

        if ((lambdaL > 0) && (lambdaR > 0)) {
          fI[nvars*p+v] = fL[nvars*p+v];
        } else if ((lambdaL < 0) && (lambdaR < 0)) {
          fI[nvars*p+v] = fR[nvars*p+v];
        } else {
          // Local Lax-Friedrichs flux for mixed/transonic case
          double alpha = max(fabs(lambdaL), fabs(lambdaR));
          fI[nvars*p+v] = 0.5 * (fL[nvars*p+v] + fR[nvars*p+v]
                               - alpha * (uR[nvars*p+v] - uL[nvars*p+v]));
        }
        */

        /* [REPLACE] Your upwinding scheme here */
        /* Default: simple averaging (non-upwinded - only for testing!) */
        fI[nvars*p+v] = 0.5 * (fL[nvars*p+v] + fR[nvars*p+v]);

      }
    }

    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
