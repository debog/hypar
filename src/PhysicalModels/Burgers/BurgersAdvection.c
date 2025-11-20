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
int BurgersAdvection( double  *a_f,   /*!< Array to hold the computed flux (same size and layout as a_u) */
                      double  *a_u,   /*!< Array containing the conserved solution */
                      int     a_dir,  /*!< Spatial dimension */
                      void    *a_s,   /*!< Solver object of type #HyPar */
                      double  a_t   /*!< Current time */
                    )
{
  HyPar     *solver = (HyPar*)   a_s;
  Burgers   *param  = (Burgers*) solver->m_physics;
  int        i, v;

  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;
  int ndims   = solver->m_ndims;
  int nvars   = solver->m_nvars;
  int index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    for (v = 0; v < nvars; v++) a_f[nvars*p+v] = 0.5 * a_u[nvars*p+v] * a_u[nvars*p+v];
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
