#include <physicalmodels/fpdoublewell.h>
#include <mpivars.h>
#include <hypar.h>

double FPDoubleWellComputeDiffNumber(void *a_s,void *a_m,double a_dt,double a_t)
{
  HyPar         *solver = (HyPar*)        a_s;
  FPDoubleWell  *params = (FPDoubleWell*) solver->m_physics;
  int           d, i, v;

  int     ndims  = solver->m_ndims;
  int     nvars  = solver->m_nvars;
  int ghosts = solver->m_ghosts;
  int     *dim   = solver->m_dim_local;
  double  *dxinv = solver->m_dxinv;

  int     offset  = 0;
  double  max_diffno = 0;
  for (d = 0; d < ndims; d++) {
    for (i = 0; i < dim[d]; i++) {
      for (v = 0; v < nvars; v++) {
        double local_diffno =   0.5 * params->m_q * a_dt
                              * dxinv[offset+ghosts+i]
                              * dxinv[offset+ghosts+i];
        if (local_diffno > max_diffno) max_diffno = local_diffno;
      }
    }
    offset += dim[d] + 2*ghosts;
  }

  return(max_diffno);
}
