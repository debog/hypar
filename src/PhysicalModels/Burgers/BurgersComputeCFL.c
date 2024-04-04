/*! @file BurgersComputeCFL.c
    @author Debojyoti Ghosh
    @brief Contains the function to compute maximum CFL over the domain for the Burgers equation.
*/

#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/burgers.h>
#include <mpivars.h>
#include <hypar.h>

/*! Computes the maximum CFL number over the domain. Note that the CFL
    is computed over the local domain on this processor only.
*/
double BurgersComputeCFL( void    *s, /*!< Solver object of type #HyPar */
                          void    *m, /*!< MPI object of type #MPIVariables */
                          double  dt, /*!< Time step size for which to compute the CFL */
                          double  t   /*!< Time */
                        )
{
  HyPar    *solver = (HyPar*)   s;
  Burgers  *params = (Burgers*) solver->physics;

  int     ndims  = solver->ndims;
  int     nvars  = solver->nvars;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;
  double  *u     = solver->u;

  int index[ndims], dir, v;

  double  max_cfl = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    for (v=0; v<nvars; v++) {
      for (dir=0; dir<ndims; dir++) {
        double dxinv;
        _GetCoordinate_(dir,index[dir],dim,ghosts,solver->dxinv,dxinv); /* 1/dx */
        double local_cfl = u[nvars*p+v]*dt*dxinv;
        if (local_cfl > max_cfl) max_cfl = local_cfl;
      }
    }
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
