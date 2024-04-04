/*! @file LinearADRComputeCFL.c
    @author Debojyoti Ghosh
    @brief Compute the maximum CFL
*/

#include <float.h>
#include <arrayfunctions.h>
#include <physicalmodels/linearadr.h>
#include <mpivars.h>
#include <hypar.h>

/*! Computes the maximum CFL number over the domain. Note that the CFL
    is computed over the local domain on this processor only.
*/
double LinearADRComputeCFL( void    *s, /*!< Solver object of type #HyPar */
                            void    *m, /*!< MPI object of type #MPIVariables */
                            double  dt, /*!< Time step size for which to compute the CFL */
                            double  t   /*!< Time */
                          )
{
  HyPar         *solver = (HyPar*)        s;
  LinearADR     *params = (LinearADR*)    solver->physics;

  int     ndims  = solver->ndims;
  int     nvars  = solver->nvars;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;

  double  max_cfl = 0;
  if (params->constant_advection == 1) {

    int d, i, v;
    for (d = 0; d < ndims; d++) {
      for (i = 0; i < dim[d]; i++) {
        for (v = 0; v < nvars; v++) {
          double dxinv; _GetCoordinate_(d,i,dim,ghosts,solver->dxinv,dxinv);
          double local_cfl = params->a[nvars*d+v]*dt*dxinv;
          if (local_cfl > max_cfl) max_cfl = local_cfl;
        }
      }
    }

  } else if (params->constant_advection == 0) {

    int d;
    for (d = 0; d < ndims; d++) {
      int index[ndims];
      int done = 0; _ArraySetValue_(index,ndims,0);
      while (!done) {
        int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
        double dxinv; _GetCoordinate_(0,index[0],dim,ghosts,solver->dxinv,dxinv);
        int v;
        for (v = 0; v < nvars; v++) {
          double a = params->a[nvars*ndims*p+nvars*d+v];
          double local_cfl = a*dt*dxinv;
          if (local_cfl > max_cfl) max_cfl = local_cfl;
        }
        _ArrayIncrementIndex_(ndims,dim,index,done);
      }
    }

  }

  return(max_cfl);
}
