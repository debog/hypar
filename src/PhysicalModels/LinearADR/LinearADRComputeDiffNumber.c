/*! @file LinearADRComputeDiffNumber.c
    @author Debojyoti Ghosh
    @brief Compute the diffusion number
*/

#include <basic.h>
#include <physicalmodels/linearadr.h>
#include <mpivars.h>
#include <hypar.h>

/*! Computes the maximum diffusion number over the domain. Note that the
 * diffusion number is computed over the local domain on this processor only.
*/
double LinearADRComputeDiffNumber( void   *s, /*!< Solver object of type #HyPar */
                                   void   *m, /*!< MPI object of type #MPIVariables */
                                   double dt, /*!< Time step size for which to compute the CFL */
                                   double t   /*!< Time */
                                 )
{
  HyPar         *solver = (HyPar*)        s;
  LinearADR     *params = (LinearADR*)    solver->physics;
  int           d, i, v;

  int     ndims  = solver->ndims;
  int     nvars  = solver->nvars;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;

  double  max_diffno = 0;
  for (d = 0; d < ndims; d++) {
    for (i = 0; i < dim[d]; i++) {
      for (v = 0; v < nvars; v++) {
        double dxinv;  _GetCoordinate_(d,i,dim,ghosts,solver->dxinv,dxinv);
        double local_diffno =   params->d[nvars*d+v] * dt * dxinv * dxinv;
        if (local_diffno > max_diffno) max_diffno = local_diffno;
      }
    }
  }

  return(max_diffno);
}
