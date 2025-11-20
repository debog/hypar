#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <physicalmodels/fpdoublewell.h>
#include <hypar.h>

int FPDoubleWellPostStep(double *a_u,void* a_s,void *a_m,double a_t,int a_iter)
{
  HyPar         *solver = (HyPar*)        a_s;
  MPIVariables  *mpi    = (MPIVariables*) a_m;
  FPDoubleWell  *params = (FPDoubleWell*) solver->m_physics;
  _DECLARE_IERR_;

  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;
  int ndims   = solver->m_ndims;
  int index[ndims];

  double local_sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double dx = 1.0 / solver->m_dxinv[index[0]+ghosts];
    local_sum     += (a_u[p] * dx);
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  double local_integral = local_sum;
  double global_integral = 0;
  IERR MPISum_double(&global_integral,&local_integral,1,&mpi->m_world); CHECKERR(ierr);
  params->m_pdf_integral = global_integral;

  return(0);
}
