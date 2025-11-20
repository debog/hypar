#include <stdlib.h>
#include <basic.h>
#include <mathfunctions.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem.h>
#include <mpivars.h>
#include <hypar.h>

double FPPowerSystemDriftFunction(int,void*,double,double,double);

double FPPowerSystemComputeCFL(void *a_s,void *a_m,double a_dt,double a_t)
{
  HyPar         *solver = (HyPar*)        a_s;
  FPPowerSystem *params = (FPPowerSystem*)solver->m_physics;

  int     ndims  = solver->m_ndims;
  int ghosts = solver->m_ghosts;
  int     *dim   = solver->m_dim_local;

  double  max_cfl = 0;
  int     index[ndims];
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    double x;     _GetCoordinate_(0,index[0],dim,ghosts,solver->m_x,x);
    double y;     _GetCoordinate_(1,index[1],dim,ghosts,solver->m_x,y);
    double dxinv; _GetCoordinate_(0,index[0],dim,ghosts,solver->m_dxinv,dxinv);
    double dyinv; _GetCoordinate_(1,index[1],dim,ghosts,solver->m_dxinv,dyinv);
    double drift_x= FPPowerSystemDriftFunction(0,params,x,y,a_t);
    double drift_y= FPPowerSystemDriftFunction(1,params,x,y,a_t);

    double local_cfl_x = absolute(drift_x) * a_dt * dxinv;
    double local_cfl_y = absolute(drift_y) * a_dt * dyinv;

    if (local_cfl_x > max_cfl) max_cfl = local_cfl_x;
    if (local_cfl_y > max_cfl) max_cfl = local_cfl_y;

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
