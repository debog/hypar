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
int LinearADRDiffusionG(  double  *f, /*!< Array to hold the computed diffusion term (same size and layout as u) */
                          double  *u, /*!< Array containing the solution */
                          int     dir,/*!< Spatial dimension (unused since this is a 1D system) */
                          void    *s, /*!< Solver object of type #HyPar */
                          double  t   /*!< Current time */
                       )
{
  HyPar     *solver = (HyPar*)     s;
  LinearADR *param  = (LinearADR*) solver->physics;
  int       i, v;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    for (v = 0; v < nvars; v++) f[nvars*p+v] = param->d[nvars*dir+v] * u[nvars*p+v];
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
int LinearADRDiffusionH(  double  *f,   /*!< Array to hold the computed diffusion term (same size and layout as u) */
                          double  *u,   /*!< Array containing the solution */
                          int     dir1, /*!< First spatial dimension of the double derivative \f$d_1\f$ */
                          int     dir2, /*!< Second spatial dimension of the double derivative \f$d_2\f$ */
                          void    *s,   /*!< Solver object of type #HyPar */
                          double  t     /*!< Current time */
                        )
{
  HyPar     *solver = (HyPar*)     s;
  LinearADR *param  = (LinearADR*) solver->physics;
  int       i, v;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int index[ndims], bounds[ndims], offset[ndims];

  if (dir1 == dir2) {

    int dir = dir1;

    /* set bounds for array index to include ghost points */
    _ArrayCopy1D_(dim,bounds,ndims);
    for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

    /* set offset such that index is compatible with ghost point arrangement */
    _ArraySetValue_(offset,ndims,-ghosts);

    int done = 0; _ArraySetValue_(index,ndims,0);
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
      for (v = 0; v < nvars; v++) f[nvars*p+v] = param->d[nvars*dir+v] * u[nvars*p+v];
      _ArrayIncrementIndex_(ndims,bounds,index,done);
    }

  } else _ArraySetValue_(f,solver->npoints_local_wghosts*nvars,0.0);


  return(0);
}
