/*! @file VlasovAdvectionCoeff.c
    @author John Loffeld
    @brief Returns the advection coefficient at a given grid index
*/

#include <float.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/vlasov.h>
#include <mpivars.h>
#include <hypar.h>

/*! Returns the advection coefficient at a given grid index
    \f$c\f$,
    where
     \f{equation}{
      c = v_i,
    \f}
    if dir < #Vlasov::ndims_x (\f$i = {\rm dir}\f$), and
     \f{equation}{
      c = E_i,
    \f}
    the electric field if dir >= #Vlasov::ndims_x (\f$i = \f$ dir-#Vlasov::ndims_x).

    \b Note: this function assumes that the electric field has already been set.
*/
double VlasovAdvectionCoeff(int*  idx, /*!< grid index */
                            int   dir, /*!< Spatial dimension */
                            void* s    /*!< Solver object of type #HyPar */
                           )
{

  HyPar  *solver = (HyPar*)  s;
  Vlasov *param  = (Vlasov*) solver->m_physics;

  int* dim    = solver->m_dim_local;
  int  ghosts = solver->m_ghosts;

  double retval = DBL_MAX;

  if (dir < param->m_ndims_x) {

    int veldim = dir + param->m_ndims_x;
    _GetCoordinate_(veldim,idx[veldim],dim,ghosts,solver->m_x,retval);

  }  else {

    int ndims_x = param->m_ndims_x;
    int dim_x[ndims_x]; _ArrayCopy1D_(dim, dim_x, ndims_x);

    int idx_x[ndims_x];
    _ArrayCopy1D_(idx, idx_x, ndims_x);
    int p; _ArrayIndex1D_(ndims_x, dim_x, idx_x, ghosts, p);
    retval = param->m_e_field[ndims_x*p+(dir-ndims_x)];

  }

  return retval;

}
