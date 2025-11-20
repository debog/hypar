/*! @file FPPowerSystem3BusDiffusion.c
    @author Debojyoti Ghosh
    @brief Compute the dissipative term for the #FPPowerSystem3Bus system.
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <hypar.h>

int FPPowerSystem3BusDissipationFunction(int,int,void*,double,double*);

/*! Compute the dissipation term for the #FPPowerSystem3Bus system */
int FPPowerSystem3BusDiffusion(
                                double  *a_f,   /*!< Array to hold the computed dissipation term vector (same layout as a_u) */
                                double  *a_u,   /*!< Array with the solution vector */
                                int     a_dir1, /*!< First spatial dimension for the dissipation term being computed */
                                int     a_dir2, /*!< Second spatial dimension for the dissipation term being computed */
                                void    *a_s,   /*!< Solver object of type #HyPar */
                                double  a_t   /*!< Current simulation time */
                              )
{
  HyPar             *solver = (HyPar*)              a_s;
  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*)  solver->m_physics;
  int               i, v;

  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;
  static int index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];
  static double dissipation[_MODEL_NDIMS_*_MODEL_NDIMS_];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_);
  for (i=0; i<_MODEL_NDIMS_; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,_MODEL_NDIMS_,-ghosts);

  int done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
    FPPowerSystem3BusDissipationFunction(a_dir1,a_dir2,params,a_t,dissipation);
    for (v = 0; v < _MODEL_NVARS_; v++) a_f[_MODEL_NVARS_*p+v] = dissipation[a_dir1*_MODEL_NDIMS_+a_dir2] * a_u[_MODEL_NVARS_*p+v];
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }

  return(0);
}
