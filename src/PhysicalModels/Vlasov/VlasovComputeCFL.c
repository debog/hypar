/*! @file VlasovComputeCFL.c
    @author John Loffeld
    @brief Contains the function to compute maximum CFL over the domain for the Vlasov equations.
*/

#include <float.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/vlasov.h>
#include <mpivars.h>
#include <hypar.h>

double VlasovAdvectionCoeff(int*, int, void*);

/*! Computes the maximum CFL number over the domain. Note that the CFL
    is computed over the local domain on this processor only.
*/
double VlasovComputeCFL( void    *s, /*!< Solver object of type #HyPar */
                         void    *m, /*!< MPI object of type #MPIVariables */
                         double  dt, /*!< Time step size for which to compute the CFL */
                         double  t   /*!< Time */
                       )
{
  HyPar   *solver = (HyPar*)  s;
  Vlasov  *params = (Vlasov*) solver->physics;

  int     ndims  = solver->ndims;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;
  double  *u     = solver->u;

  double max_cfl = 0;
  int done = 0;
  int index[ndims];
  _ArraySetValue_(index,ndims,0);

  while (!done) {

    for (int dir=0; dir<ndims; dir++) {
      double dxinv;
      _GetCoordinate_(dir,index[dir],dim,ghosts,solver->dxinv,dxinv);
      double eig = VlasovAdvectionCoeff(index, dir, solver);
      double local_cfl = eig*dt*dxinv;
      if (local_cfl > max_cfl) max_cfl = local_cfl;
    }

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
