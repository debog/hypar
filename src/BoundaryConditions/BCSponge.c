/*! @file BCSponge.c
    @author Debojyoti Ghosh
    @brief Sponge boundary condition
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

/*! Applies the sponge boundary condition: This function computes the source term required to apply a
    sponge boundary condition that gradually relaxes the solution to a specified state. This boundary
    condition is different from other boundary conditions in the sense that it is applied at interior
    grid points (but within the defined sponge zone).
    \n\n
    The source term for the sponge is computed as:
    \f{align}{
      {\bf S}_i &= \sigma_i \left( {\bf U}_i - {\bf U}_{\rm ref} \right),\\
      \sigma_i &= \frac {x_i - x_{\rm start}} {x_{\rm end} - x_{\rm start}}
    \f}
    where \f$i\f$ is the grid index along the spatial dimension of the sponge, \f${\bf U}_{\rm ref}\f$
    is the specified state to which the solution is relaxed, and \f$x_i\f$, \f$x_{\rm start}\f$, and
    \f$x_{\rm end}\f$ are the spatial coordinates of the grid point, sponge start, and sponge end,
    respectively, along the spatial dimension of the sponge.
*/
int BCSpongeSource(
                    void    *a_b,     /*!< Boundary object of type #DomainBoundary */
                    int     a_ndims,  /*!< Number of spatial dimensions */
                    int     a_nvars,  /*!< Number of variables/DoFs per grid point */
                    int     a_ghosts, /*!< Number of ghost points */
                    int     *a_size,  /*!< Integer array with the number of grid points in each spatial dimension */
                    double  *a_grid,  /*!< 1D array with the spatial coordinates of the grid points, one dimension after the other */
                    double  *a_u,     /*!< Solution */
                    double  *a_source /*!< Source term to which the sponge term is added */
                  )
{
  DomainBoundary *boundary = (DomainBoundary*) a_b;
  int            dim       = boundary->m_dim;
  int            face      = boundary->m_face;
  double         *uref     = boundary->m_SpongeValue;
  double         *xmin     = boundary->m_xmin;
  double         *xmax     = boundary->m_xmax;
  int            v;

  if (boundary->m_on_this_proc) {
    int bounds[a_ndims], indexb[a_ndims];
    _ArraySubtract1D_(bounds,boundary->m_ie,boundary->m_is,a_ndims);
    _ArraySetValue_(indexb,a_ndims,0);
    int done = 0;
    while (!done) {
      int i = indexb[dim] + boundary->m_is[dim];
      double x, xstart, xend;
      _GetCoordinate_(dim,i,a_size,a_ghosts,a_grid,x);
      xstart = xmin[dim];
      xend   = xmax[dim];
      /* calculate sigma */
      double sigma;
      if (face > 0) sigma = (x - xstart) / (xend - xstart);
      else          sigma = (x - xend  ) / (xstart - xend);
      /* add to the a_source term */
      int p; _ArrayIndex1DWO_(a_ndims,a_size,indexb,boundary->m_is,a_ghosts,p);
      for (v=0; v<a_nvars; v++) a_source[a_nvars*p+v] -= (sigma * (a_u[a_nvars*p+v]-uref[v]));
      _ArrayIncrementIndex_(a_ndims,bounds,indexb,done);
    }
  }
  return(0);
}

/*! Dummy function to ensure consistency with the overall boundary condition
    implementation. The actual sponge boundary condition is implemented by
    BCSpongeSource()
*/
int BCSpongeUDummy(
                    void    *a_b,     /*!< Boundary object of type #DomainBoundary */
                    void    *a_m,     /*!< MPI object of type #MPIVariables */
                    int     a_ndims,  /*!< Number of spatial dimensions */
                    int     a_nvars,  /*!< Number of variables/DoFs per grid point */
                    int     *a_size,  /*!< Integer array with the number of grid points in each spatial dimension */
                    int     a_ghosts, /*!< Number of ghost points */
                    double  *a_phi,   /*!< The solution array on which to apply the boundary condition */
                    double  a_waqt    /*!< Current solution time */
                  )
{
  return(0);
}
