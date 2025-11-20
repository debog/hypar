/*! @file BCReflect.c
    @author Debojyoti Ghosh
    @brief Reflection boundary condition
*/
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

/*! Applies the reflection boundary condition: The values at the physical
    boundary ghost points are set to the negative of the interior values
    adjacent to the boundary.
*/
int BCReflectU(
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
  DomainBoundary *boundary = (DomainBoundary*) a_b;

  int dim   = boundary->m_dim;
  int face  = boundary->m_face;

  if (boundary->m_on_this_proc) {
    int bounds[a_ndims], indexb[a_ndims], indexi[a_ndims];
    _ArraySubtract1D_(bounds,boundary->m_ie,boundary->m_is,a_ndims);
    _ArraySetValue_(indexb,a_ndims,0);
    int done = 0;
    while (!done) {
      int p1, p2;
      _ArrayCopy1D_(indexb,indexi,a_ndims);
      _ArrayAdd1D_(indexi,indexi,boundary->m_is,a_ndims);
      if      (face ==  1) indexi[dim] = a_ghosts-1-indexb[dim];
      else if (face == -1) indexi[dim] = a_size[dim]-indexb[dim]-1;
      else return(1);
      _ArrayIndex1DWO_(a_ndims,a_size,indexb,boundary->m_is,a_ghosts,p1);
      _ArrayIndex1D_(a_ndims,a_size,indexi,a_ghosts,p2);
      _ArrayScaleCopy1D_((a_phi+a_nvars*p2),(-1.0),(a_phi+a_nvars*p1),a_nvars);
      _ArrayIncrementIndex_(a_ndims,bounds,indexb,done);
    }
  }
  return(0);
}
