/*! @file BCPeriodic.c
    @author Debojyoti Ghosh
    @brief Periodic boundary conditions
*/
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <boundaryconditions.h>

/*! Applies periodic boundary conditions: Implemented by copying the solution
    from the other end of the domain into the physical boundary ghost points.
    \n\n
    **Note**: This function only acts if the the number of processors is 1 along
    the spatial dimension this boundary corresponds to. If there are more than 1
    processors along this dimension, periodicity is handled by MPIExchangeBoundariesnD()
    to minimize communication.
*/
int BCPeriodicU(
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
  MPIVariables   *mpi      = (MPIVariables*)   a_m;

  int dim   = boundary->m_dim;
  int face  = boundary->m_face;

  if ((boundary->m_on_this_proc) && (mpi->m_iproc[dim] == 1)) {
    int bounds[a_ndims], index1[a_ndims], index2[a_ndims];
    _ArraySubtract1D_(bounds,boundary->m_ie,boundary->m_is,a_ndims);
    _ArraySetValue_(index1,a_ndims,0);
    _ArraySetValue_(index2,a_ndims,0);
    int done = 0;
    while (!done) {
      int p1 = 0, p2 = 0;
      _ArrayCopy1D_(index1,index2,a_ndims);
      if (face == 1) {
        index2[dim] = index1[dim] + a_size[dim]-a_ghosts;
        _ArrayIndex1DWO_(a_ndims,a_size,index1,boundary->m_is,a_ghosts,p1);
        _ArrayIndex1D_(a_ndims,a_size,index2,a_ghosts,p2);
      } else if (face == -1) {
        _ArrayIndex1DWO_(a_ndims,a_size,index1,boundary->m_is,a_ghosts,p1);
        _ArrayIndex1D_(a_ndims,a_size,index1,a_ghosts,p2);
      }
      _ArrayCopy1D_((a_phi+a_nvars*p2),(a_phi+a_nvars*p1),a_nvars);
      _ArrayIncrementIndex_(a_ndims,bounds,index1,done);
    }
  }

  return(0);
}
