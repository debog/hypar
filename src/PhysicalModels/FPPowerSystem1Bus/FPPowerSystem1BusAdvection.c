#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem1bus.h>
#include <hypar.h>

int FPPowerSystem1BusAdvection(double *a_f,double *a_u,int a_dir,void *a_s,double a_t)
{
  HyPar *solver = (HyPar*) a_s;
  int   *dim    = solver->m_dim_local;
  int   ghosts  = solver->m_ghosts;
  int   ndims   = solver->m_ndims;
  int   nvars   = solver->m_nvars;

  /* calculate total size of arrays */
  int bounds[ndims]; _ArrayAddCopy1D_(dim,(2*ghosts),bounds,ndims);
  int size;          _ArrayProduct1D_(bounds,ndims,size); size *= nvars;

  /* a_f = a_u */
  _ArrayCopy1D_(a_u,a_f,size);

  return(0);
}
