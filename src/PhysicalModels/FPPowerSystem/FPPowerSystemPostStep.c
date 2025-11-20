#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <physicalmodels/fppowersystem.h>
#include <hypar.h>

int FPPowerSystemPostStep(double *a_u,void* a_s,void *a_m,double a_t,int a_iter)
{
  HyPar         *solver = (HyPar*)         a_s;
  MPIVariables  *mpi    = (MPIVariables*)  a_m;
  FPPowerSystem *params = (FPPowerSystem*) solver->m_physics;
  _DECLARE_IERR_;

  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;
  int ndims   = solver->m_ndims;
  int index[ndims];

  double local_sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    double dxinv; _GetCoordinate_(0,index[0],dim,ghosts,solver->m_dxinv,dxinv);
    double dyinv; _GetCoordinate_(1,index[1],dim,ghosts,solver->m_dxinv,dyinv);
    local_sum     += (a_u[p] / (dxinv * dyinv));
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  double local_integral = local_sum;
  double global_integral = 0;
  IERR MPISum_double(&global_integral,&local_integral,1,&mpi->m_world); CHECKERR(ierr);
  params->m_pdf_integral = global_integral;

  return(0);
}
