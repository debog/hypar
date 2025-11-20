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
double LinearADRComputeCFL( void    *a_s, /*!< Solver object of type #HyPar */
                            void    *a_m, /*!< MPI object of type #MPIVariables */
                            double  a_dt, /*!< Time step size for which to compute the CFL */
                            double  a_t   /*!< Time */
                          )
{
  HyPar         *solver = (HyPar*)        a_s;
  LinearADR     *params = (LinearADR*)    solver->m_physics;

  int     ndims  = solver->m_ndims;
  int     nvars  = solver->m_nvars;
  int ghosts = solver->m_ghosts;
  int     *dim   = solver->m_dim_local;

  double  max_cfl = 0;
  if (params->m_constant_advection == 1) {

    int d, i, v;
    for (d = 0; d < ndims; d++) {
      for (i = 0; i < dim[d]; i++) {
        for (v = 0; v < nvars; v++) {
          double dxinv; _GetCoordinate_(d,i,dim,ghosts,solver->m_dxinv,dxinv);
          double local_cfl = params->m_a[nvars*d+v]*a_dt*dxinv;
          if (local_cfl > max_cfl) max_cfl = local_cfl;
        }
      }
    }

  } else if (params->m_constant_advection == 0) {

    int d;
    for (d = 0; d < ndims; d++) {
      int index[ndims];
      int done = 0; _ArraySetValue_(index,ndims,0);
      while (!done) {
        int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
        double dxinv; _GetCoordinate_(0,index[0],dim,ghosts,solver->m_dxinv,dxinv);
        int v;
        for (v = 0; v < nvars; v++) {
          double a = params->m_a[nvars*ndims*p+nvars*d+v];
          double local_cfl = a*a_dt*dxinv;
          if (local_cfl > max_cfl) max_cfl = local_cfl;
        }
        _ArrayIncrementIndex_(ndims,dim,index,done);
      }
    }

  }

  return(max_cfl);
}
