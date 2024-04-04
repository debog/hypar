/*! @file LinearADRAdvection.c
    @author Debojyoti Ghosh
    @brief Function to evaluate the advection term in the
           linear advection-diffusion-reaction model
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/linearadr.h>
#include <hypar.h>

/*! Evaluate the advection term in the linear advection-diffusion-reaction model:\n

    Compute
    \f{equation}
    a_d u
    \f}
    given \f$u\f$ and \f$d\f$ in the hyperbolic term
    \f{equation}{
    \sum_d \frac {\partial} {\partial x_d} \left( a_d u \right)
    \f}

    + If the advection field is spatially varying, then the advection speed
      array (LinearADR::a) is expected to have the size (LinearADR::adv_arr_size):
      (number of grid points with ghosts - HyPar::npoints_local_wghosts)
      X (number of spatial dimensions - HyPar::ndims)
      X (number of solution components - HyPar::nvars),
      with the innermost loop being the solution components and outermost loop
      being the grid points.
    + If the advection field is constant, then the advection speed
      array (LinearADR::a) is expected to have the size (LinearADR::adv_arr_size):
      (number of spatial dimensions - HyPar::ndims)
      X (number of solution components - HyPar::nvars),
      with the inner loop being the solution components.
*/
int LinearADRAdvection( double  *f, /*!< Array to hold the computed flux (same size and layout as u) */
                        double  *u, /*!< Array containing the solution */
                        int     dir,/*!< Spatial dimension \f$d\f$ */
                        void    *s, /*!< Solver object of type #HyPar */
                        double  t   /*!< Current time */
                      )
{
  HyPar     *solver = (HyPar*)     s;
  LinearADR *param  = (LinearADR*) solver->physics;
  double    *adv    = param->a;
  int       i, v;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars;

  int index[ndims],
      bounds[ndims],
      offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  if (param->constant_advection == 1) {
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
      for (v = 0; v < nvars; v++) {
        f[nvars*p+v] = param->a[nvars*dir+v] * u[nvars*p+v];
      }
      _ArrayIncrementIndex_(ndims,bounds,index,done);
    }
  } else if (param->constant_advection == 0) {
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
      for (v = 0; v < nvars; v++) {
        double a = adv[nvars*ndims*p+nvars*dir+v];
        f[nvars*p+v] = a * u[nvars*p+v];
      }
      _ArrayIncrementIndex_(ndims,bounds,index,done);
    }
  } else {
    while (!done) {
      int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
      for (v = 0; v < nvars; v++) {
        f[nvars*p+v] = 0.0;
      }
      _ArrayIncrementIndex_(ndims,bounds,index,done);
    }
  }

  return(0);
}
