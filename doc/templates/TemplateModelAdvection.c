/*! @file TemplateModelAdvection.c
    @author [YOUR NAME]
    @brief Compute the hyperbolic flux for the Template Model
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/template_model.h>
#include <hypar.h>

/*!
 * Compute the hyperbolic flux function over the local domain.
 *
 * [REPLACE] Mathematical description:
 * \f{equation}{
 *   {\bf f}({\bf u}) = [REPLACE with your flux function]
 * \f}
 *
 * This function computes the flux at each grid point (including ghost points).
 * The flux array \a f has the same size and layout as the solution array \a u.
*/
int TemplateModelAdvection(
  double  *f,     /*!< Array to hold the computed flux (same size and layout as u) */
  double  *u,     /*!< Array containing the conserved solution */
  int     dir,    /*!< Spatial dimension (0 for x, 1 for y, 2 for z) */
  void    *s,     /*!< Solver object of type #HyPar */
  double  t       /*!< Current simulation time */
)
{
  HyPar         *solver  = (HyPar*)         s;
  TemplateModel *physics = (TemplateModel*) solver->m_physics;

  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;
  int ndims   = solver->m_ndims;
  int nvars   = solver->m_nvars;

  /*
   * Array indexing in HyPar:
   * - Solution u is a 1D array: u[nvars*p + v]
   *   where p is the point index and v is the variable index
   * - Use macros for multi-dimensional to 1D index conversion
  */

  int index[ndims], bounds[ndims], offset[ndims];
  int i, v;

  /* Set bounds to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* Set offset for ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  /* Loop over all grid points including ghost points */
  int done = 0;
  _ArraySetValue_(index,ndims,0);
  while (!done) {

    /* Calculate 1D index from multi-dimensional index */
    int p;
    _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);

    /* [REPLACE] Compute flux for each variable at this point */
    /* Example for scalar equation: f = a * u */
    /*
    for (v = 0; v < nvars; v++) {
      f[nvars*p+v] = physics->param1 * u[nvars*p+v];
    }
    */

    /* [REPLACE] Example for system of equations (e.g., Euler) */
    /*
    double rho  = u[nvars*p+0];
    double rhou = u[nvars*p+1];
    double E    = u[nvars*p+2];
    double u_vel = rhou/rho;
    double p = (physics->param1 - 1.0) * (E - 0.5*rho*u_vel*u_vel);

    if (dir == 0) {  // x-direction flux
      f[nvars*p+0] = rhou;
      f[nvars*p+1] = rhou*u_vel + p;
      f[nvars*p+2] = (E + p) * u_vel;
    }
    */

    /* [REPLACE] Your flux calculation here */
    for (v = 0; v < nvars; v++) {
      f[nvars*p+v] = 0.0;  // PLACEHOLDER - implement your flux
    }

    /* Increment to next grid point */
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
