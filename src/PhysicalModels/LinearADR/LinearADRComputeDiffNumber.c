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
double LinearADRComputeDiffNumber( void   *a_s, /*!< Solver object of type #HyPar */
                                   void   *a_m, /*!< MPI object of type #MPIVariables */
                                   double a_dt, /*!< Time step size for which to compute the CFL */
                                   double  a_t   /*!< Time */
                                 )
{
  HyPar         *solver = (HyPar*)        a_s;
  LinearADR     *params = (LinearADR*)    solver->m_physics;
  int           d, i, v;

  int     ndims  = solver->m_ndims;
  int     nvars  = solver->m_nvars;
  int ghosts = solver->m_ghosts;
  int     *dim   = solver->m_dim_local;

  double  max_diffno = 0;
  for (d = 0; d < ndims; d++) {
    for (i = 0; i < dim[d]; i++) {
      for (v = 0; v < nvars; v++) {
        double dxinv;  _GetCoordinate_(d,i,dim,ghosts,solver->m_dxinv,dxinv);
        double local_diffno =   params->m_d[nvars*d+v] * a_dt * dxinv * dxinv;
        if (local_diffno > max_diffno) max_diffno = local_diffno;
      }
    }
  }

  return(max_diffno);
}
