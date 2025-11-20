#include <basic.h>
#include <mathfunctions.h>
#include <physicalmodels/fpdoublewell.h>
#include <mpivars.h>
#include <hypar.h>

double FPDoubleWellComputeCFL(void *a_s,void *a_m,double a_dt,double a_t)
{
  HyPar         *solver = (HyPar*)        a_s;
  int           d, i, v;

  int     ndims  = solver->m_ndims;
  int     nvars  = solver->m_nvars;
  int ghosts = solver->m_ghosts;
  int     *dim   = solver->m_dim_local;
  double  *dxinv = solver->m_dxinv;

  int     offset  = 0;
  double  max_cfl = 0;
  for (d = 0; d < ndims; d++) {
    for (i = 0; i < dim[d]; i++) {
      for (v = 0; v < nvars; v++) {
        double x; _GetCoordinate_(0,i,dim,ghosts,solver->m_x,x);
        double local_cfl =  absolute(drift(x)) * a_dt
                          * dxinv[offset+ghosts+i];
        if (local_cfl > max_cfl) max_cfl = local_cfl;
      }
    }
    offset += dim[d] + 2*ghosts;
  }

  return(max_cfl);
}
