/*! @file TemplateModelComputeCFL.c
    @author [YOUR NAME]
    @brief Compute the maximum CFL number for the Template Model
*/

#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/template_model.h>
#include <mpivars.h>
#include <hypar.h>

/*! 
 * Compute the maximum CFL number over the local domain.
 * 
 * The CFL number is defined as:
 * \f{equation}{
 *   CFL = \lambda \frac{\Delta t}{\Delta x}
 * \f}
 * where \f$\lambda\f$ is the maximum wave speed (eigenvalue).
 * 
 * For explicit time-stepping stability, CFL should typically be < 1.0
 * (exact limit depends on the time integration scheme).
 * 
 * Note: This function computes CFL only over the local domain on this
 * processor. The global maximum across all processors is computed separately
 * by the main solver.
*/
double TemplateModelComputeCFL(
  void    *s,     /*!< Solver object of type #HyPar */
  void    *m,     /*!< MPI object of type #MPIVariables */
  double  dt,     /*!< Time step size for which to compute the CFL */
  double  t       /*!< Current simulation time */
)
{
  HyPar         *solver  = (HyPar*)         s;
  TemplateModel *physics = (TemplateModel*) solver->physics;

  int     ndims  = solver->ndims;
  int     nvars  = solver->nvars;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;
  double  *u     = solver->u;

  double max_cfl = 0.0;
  int index[ndims];
  int done = 0;

  _ArraySetValue_(index,ndims,0);

  /* Loop over all interior grid points (not ghost points) */
  while (!done) {
    
    /* Calculate 1D index */
    int p;
    _ArrayIndex1D_(ndims,dim,index,ghosts,p);

    /* [REPLACE] Compute maximum wave speed at this point */
    /* For each variable and direction, compute local CFL */
    
    int v, dir;
    for (v=0; v<nvars; v++) {
      for (dir=0; dir<ndims; dir++) {
        
        /* Get inverse grid spacing in this direction */
        double dxinv;
        _GetCoordinate_(dir,index[dir],dim,ghosts,solver->dxinv,dxinv);

        /* [REPLACE] Calculate wave speed (eigenvalue) */
        
        /* Example 1: Constant advection speed */
        /*
        double wave_speed = physics->param1;
        */

        /* Example 2: Variable wave speed depending on solution */
        /*
        double wave_speed = u[nvars*p+v];
        */

        /* Example 3: For systems, compute maximum eigenvalue */
        /*
        double rho = u[nvars*p+0];
        double vel = u[nvars*p+1] / rho;
        double c = sqrt(physics->param1 * physics->param2 / rho);  // sound speed
        double wave_speed = fabs(vel) + c;  // maximum characteristic speed
        */

        /* [REPLACE] Your wave speed calculation */
        double wave_speed = 1.0;  // PLACEHOLDER

        /* Calculate local CFL number */
        double local_cfl = fabs(wave_speed) * dt * dxinv;

        /* Update maximum */
        if (local_cfl > max_cfl) max_cfl = local_cfl;
      }
    }

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
