/*! @file VlasovComputeCFL.c
    @author Debojyoti Ghosh
    @brief Contains the function to compute maximum CFL over the domain for the Vlasov equations.
*/

#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/vlasov.h>
#include <mpivars.h>
#include <hypar.h>

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
  int     nvars  = solver->nvars;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;
  double  *u     = solver->u;

#ifdef fftw
  double *field = params->field;
#endif

  bool self_consistent_electric_field = params->self_consistent_electric_field;

  int index[ndims], dir, v;

  double  max_cfl = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {

    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);

    for (dir=0; dir<ndims; dir++) {

      double dxinv;
      _GetCoordinate_(dir,index[dir],dim,ghosts,solver->dxinv,dxinv); /* 1/dx */

      double eig;
      if (dir == 0) {
        /* Velocity coordinate is the velocity */
        _GetCoordinate_(1,index[1],dim,ghosts,solver->x,eig);
      } else {
        if (self_consistent_electric_field) {
          /* assumes field has been calculated just prior in VlasovAdvection */
          eig = field[index[0]];
        } else {
          /* Prescribed electric field is the velocity */
          double x;
          _GetCoordinate_(0,index[0],dim,ghosts,solver->x,x);
          eig = 0.1 * cos(x);
        }
      }

      double local_cfl = eig*dt*dxinv; 
      if (local_cfl > max_cfl) max_cfl = local_cfl;

    }
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
