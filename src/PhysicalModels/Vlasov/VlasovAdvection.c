/*! @file VlasovAdvection.c
    @author John Loffeld
    @brief Contains the function to compute the hyperbolic flux for the Vlasov equations over the domain.
*/

#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/vlasov.h>
#include <mpivars.h>
#include <hypar.h>

double VlasovAdvectionCoeff(int*, int, void*);

/*! Compute the hyperbolic flux over the local domain.\n
    \f{equation}{
      {\bf F}\left({\bf u}(x,v)\right) = c {\bf u}(x,v),
    \f}
    where the advection coefficient \f$c\f$ is computed
    by VlasovAdvectionCoeff().
*/
int VlasovAdvection( double  *f,   /*!< Array to hold the computed flux (same size and layout as u) */
                     double  *u,   /*!< Array containing the conserved solution */
                     int      dir, /*!< Spatial dimension */
                     void    *s,   /*!< Solver object of type #HyPar */
                     double   t    /*!< Current time */
                   )
{

  HyPar  *solver = (HyPar*)  s;
  Vlasov *param  = (Vlasov*) solver->physics;

  int* dim    = solver->dim_local;
  int  ghosts = solver->ghosts;
  int  ndims  = solver->ndims;

  // set bounds for array index to include ghost points
  int bounds[ndims];
  _ArrayCopy1D_(dim,bounds,ndims);
  for (int i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  // set offset such that index is compatible with ghost point arrangement
  int offset[ndims];
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0;
  int index_wog[ndims];
  int index[ndims]; _ArraySetValue_(index,ndims,0);
  while (!done) {

    _ArrayCopy1D_(index, index_wog, ndims);
    for (int i = 0; i < ndims; i++) index_wog[i] -= ghosts;
    double adv_coeff = VlasovAdvectionCoeff(index_wog, dir, solver);

    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    f[p] = adv_coeff * u[p];
    _ArrayIncrementIndex_(ndims,bounds,index,done);

  }

  return 0;
}
