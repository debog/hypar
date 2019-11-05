/*! @file BurgersAdvection.c
    @author John Loffeld
    @brief Contains the function to compute the hyperbolic flux for the Burgers equations over the domain.
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/burgers.h>
#include <hypar.h>

/*! Compute the hyperbolic flux over the local domain.\n
    \f{equation}{
      {\bf F}\left({\bf u}\right) = 0.5 {\bf u}^2
    \f}
*/
int BurgersAdvection( double  *f,   /*!< Array to hold the computed flux (same size and layout as u) */
                      double  *u,   /*!< Array containing the conserved solution */
                      int     dir,  /*!< Spatial dimension */
                      void    *s,   /*!< Solver object of type #HyPar */
                      double  t     /*!< Current time */
                    )
{
  HyPar     *solver = (HyPar*)   s;
  Burgers   *param  = (Burgers*) solver->physics;
  int        i, v;

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
    for (v = 0; v < nvars; v++) f[nvars*p+v] = 0.5 * u[nvars*p+v] * u[nvars*p+v];
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
