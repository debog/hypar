/*! @file LinearADRDiffusion.c
    @author Debojyoti Ghosh
    @brief Function to evaluate the diffusion term in the
           linear advection-diffusion-reaction model
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/linearadr.h>
#include <mpivars.h>
#include <hypar.h>

/*! Evaluate the diffusion term in the linear advection-diffusion-reaction model
    for a "pure Laplacian" type operator (no cross derivatives):\n

    Compute
    \f{equation}
    \nu_d u
    \f}
    given \f$u\f$ and \f$d\f$ in the parabolic term
    \f{equation}{
    \sum_d \frac {\partial^2} {\partial x_d^2} \left( \nu_d u \right)
    \f}
*/
int LinearADRDiffusionG(  double  *a_f, /*!< Array to hold the computed diffusion term (same size and layout as a_u) */
                          double  *a_u, /*!< Array containing the solution */
                          int     a_dir,/*!< Spatial dimension (unused since this is a 1D system) */
                          void    *a_s, /*!< Solver object of type #HyPar */
                          double  a_t   /*!< Current time */
                       )
{
  HyPar     *solver = (HyPar*)     a_s;
  LinearADR *param  = (LinearADR*) solver->m_physics;
  int       i, v;

  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;
  int ndims   = solver->m_ndims;
  int nvars   = solver->m_nvars;
  int index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    for (v = 0; v < nvars; v++) a_f[nvars*p+v] = param->m_d[nvars*a_dir+v] * a_u[nvars*p+v];
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}

/*! Evaluate the diffusion term in the linear advection-diffusion-reaction model
    for a parabolic operator with no cross derivatives:\n

    Compute
    \f{equation}
    \nu_d u
    \f}
    given \f$u\f$ and \f$d_1,d_2\f$ in the parabolic term
    \f{equation}{
    \sum_{d_1}\sum_{d_2} \frac {\partial^2} {\partial x_{d_1} \partial x_{d_2}} \left( \nu_d u \right)
    \f}

    \b Note: it's not correctly implemented. Will implement when necessary.
*/
int LinearADRDiffusionH(  double  *a_f,   /*!< Array to hold the computed diffusion term (same size and layout as a_u) */
                          double  *a_u,   /*!< Array containing the solution */
                          int     a_dir1, /*!< First spatial dimension of the double derivative \f$d_1\f$ */
                          int     a_dir2, /*!< Second spatial dimension of the double derivative \f$d_2\f$ */
                          void    *a_s,   /*!< Solver object of type #HyPar */
                          double  a_t   /*!< Current time */
                        )
{
  HyPar     *solver = (HyPar*)     a_s;
  LinearADR *param  = (LinearADR*) solver->m_physics;
  int       i, v;

  int *dim    = solver->m_dim_local;
  int ghosts  = solver->m_ghosts;
  int ndims   = solver->m_ndims;
  int nvars   = solver->m_nvars;
  int index[ndims], bounds[ndims], offset[ndims];

  if (a_dir1 == a_dir2) {

    int dir = a_dir1;

    /* set bounds for array index to include ghost points */
    _ArrayCopy1D_(dim,bounds,ndims);
    for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

    /* set offset such that index is compatible with ghost point arrangement */
    _ArraySetValue_(offset,ndims,-ghosts);

    int done = 0; _ArraySetValue_(index,ndims,0);
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
      for (v = 0; v < nvars; v++) a_f[nvars*p+v] = param->m_d[nvars*dir+v] * a_u[nvars*p+v];
      _ArrayIncrementIndex_(ndims,bounds,index,done);
    }

  } else _ArraySetValue_(a_f,solver->m_npoints_local_wghosts*nvars,0.0);


  return(0);
}
