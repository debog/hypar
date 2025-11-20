/*! @file BCDirichlet.c
    @author Debojyoti Ghosh
    @brief Dirichlet boundary conditions
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>

/*! Applies (steady) Dirichlet boundary conditions for the solution: the ghost points at the physical
    boundaries are set to specified values */
int BCDirichletU(
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

  if (boundary->m_on_this_proc) {
    int bounds[a_ndims], indexb[a_ndims];
    _ArraySubtract1D_(bounds,boundary->m_ie,boundary->m_is,a_ndims);
    _ArraySetValue_(indexb,a_ndims,0);
    int done = 0;
    while (!done) {
      int p; _ArrayIndex1DWO_(a_ndims,a_size  ,indexb,boundary->m_is,a_ghosts,p);
      _ArrayCopy1D_((boundary->m_DirichletValue),(a_phi+a_nvars*p),a_nvars);
      _ArrayIncrementIndex_(a_ndims,bounds,indexb,done);
    }
  }
  return(0);
}
