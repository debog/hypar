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
double VlasovAdvectionCoeff(int*  a_idx, /*!< grid index */
                            int   a_dir, /*!< Spatial dimension */
                            void* a_s    /*!< Solver object of type #HyPar */
                           )
{

  HyPar  *solver = (HyPar*)  a_s;
  Vlasov *param  = (Vlasov*) solver->m_physics;

  int* dim    = solver->m_dim_local;
  int  ghosts = solver->m_ghosts;

  double retval = DBL_MAX;

  if (a_dir < param->m_ndims_x) {

    int veldim = a_dir + param->m_ndims_x;
    _GetCoordinate_(veldim,a_idx[veldim],dim,ghosts,solver->m_x,retval);

  }  else {

    int ndims_x = param->m_ndims_x;
    int dim_x[ndims_x]; _ArrayCopy1D_(dim, dim_x, ndims_x);

    int idx_x[ndims_x];
    _ArrayCopy1D_(a_idx, idx_x, ndims_x);
    int p; _ArrayIndex1D_(ndims_x, dim_x, idx_x, ghosts, p);
    retval = param->m_e_field[ndims_x*p+(a_dir-ndims_x)];

  }

  return retval;

}
