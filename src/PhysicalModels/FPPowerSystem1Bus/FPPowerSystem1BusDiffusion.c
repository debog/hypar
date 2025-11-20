#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem1bus.h>
#include <hypar.h>

double FPPowerSystem1BusDissipationFunction(int,int,void*,double);

int FPPowerSystem1BusDiffusionLaplacian(double *a_f,double *a_u,int a_dir,void *a_s,double a_t)
{
  HyPar             *solver = (HyPar*)             a_s;
  FPPowerSystem1Bus *params = (FPPowerSystem1Bus*) solver->m_physics;

  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;
  int ndims   = solver->m_ndims;
  int nvars   = solver->m_nvars;

  /* calculate total size of arrays */
  int bounds[ndims]; _ArrayAddCopy1D_(dim,(2*ghosts),bounds,ndims);
  int size;          _ArrayProduct1D_(bounds,ndims,size); size *= nvars;

  /* calculate dissipation coefficient  -- constant in x and y */
  double dissipation = FPPowerSystem1BusDissipationFunction(a_dir,a_dir,params,a_t);

  /* a_f = dissipation * a_u */
  _ArrayScaleCopy1D_(a_u,dissipation,a_f,size);

  return(0);
}

int FPPowerSystem1BusDiffusionGeneral(double *a_f,double *a_u,int a_dir1,int a_dir2,void *a_s,double a_t)
{
  HyPar             *solver = (HyPar*)             a_s;
  FPPowerSystem1Bus *params = (FPPowerSystem1Bus*) solver->m_physics;

  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;
  int ndims   = solver->m_ndims;
  int nvars   = solver->m_nvars;

  /* calculate total size of arrays */
  int bounds[ndims]; _ArrayAddCopy1D_(dim,(2*ghosts),bounds,ndims);
  int size;          _ArrayProduct1D_(bounds,ndims,size); size *= nvars;

  /* calculate dissipation coefficient  -- constant in x and y */
  double dissipation = FPPowerSystem1BusDissipationFunction(a_dir1,a_dir2,params,a_t);

  /* a_f = dissipation * a_u */
  _ArrayScaleCopy1D_(a_u,dissipation,a_f,size);

  return(0);
}
