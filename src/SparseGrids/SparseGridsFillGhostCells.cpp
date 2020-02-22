/*! @file SparseGridsFillGhostCells.cpp
    @brief Fill the ghost cells of a solution
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <arrayfunctions.h>
#include <sparse_grids_simulation.h>

/*! Fill the ghost cells of a solution. Note that the solution array *must* be
    a global one (not one that is partitioned among MPI ranks). This is not a
    parallel operation (it will execute independently on multiple MPI ranks, if
    called by multiple processes).
  
    For periodicity along any dimension, the ghost cells are filled appropriately
    from the other side of the domain. Otherwise, the interior data is extrapolated
    by a 4th order polynomial (assuming uniform grid spacing).
*/
void SparseGridsSimulation::fillGhostCells( const GridDimensions& a_dim, /*!< grid dimensions of solution */
                                            const int a_ngpt, /*!< number of ghost cells */
                                            double* const  a_u, /*!< solution array */
                                            const int a_nvars /*!< number of vector components of the solution */
                                          )
{
  if (m_is_periodic.size() != m_ndims) {
    fprintf(stderr, "Error in SparseGridsSimulation::fillGhostCells() -\n");
    fprintf(stderr, "m_is_periodic.size() != m_ndims \n");
    exit(1);
  }

  for (int d = 0; d < m_ndims; d++) {

    if (m_is_periodic[d]) {
      
      /* periodic case */

      int bounds[m_ndims];
      _ArrayCopy1D_(a_dim, bounds, m_ndims);
      bounds[d] = a_ngpt;

      int index[m_ndims];
      _ArraySetValue_(index, m_ndims, 0);

      int done = 0;
      while (!done) {

        {
          /* low end - face = 1 */

          int p_gpt = 0, 
              p_int = 0;

          int offset_gpt[m_ndims];
          _ArraySetValue_(offset_gpt, m_ndims, 0);
          offset_gpt[d] = -a_ngpt;

          int offset_int[m_ndims];
          _ArraySetValue_(offset_int, m_ndims, 0);
          offset_int[d] = a_dim[d] - a_ngpt;

          _ArrayIndex1DWO_(m_ndims, a_dim, index, offset_gpt, a_ngpt, p_gpt);
          _ArrayIndex1DWO_(m_ndims, a_dim, index, offset_int, a_ngpt, p_int);

          _ArrayCopy1D_((a_u+a_nvars*p_int), (a_u+a_nvars*p_gpt), a_nvars);
        }

        {
          /* high end - face = -1 */

          int p_gpt = 0, 
              p_int = 0;

          int offset_gpt[m_ndims];
          _ArraySetValue_(offset_gpt, m_ndims, 0);
          offset_gpt[d] = a_dim[d];

          _ArrayIndex1DWO_(m_ndims, a_dim, index, offset_gpt, a_ngpt, p_gpt);
          _ArrayIndex1D_(m_ndims, a_dim, index, a_ngpt, p_int);

          _ArrayCopy1D_((a_u+a_nvars*p_int), (a_u+a_nvars*p_gpt), a_nvars);
        }

        _ArrayIncrementIndex_(m_ndims, bounds, index, done);

      }

    } else {

      /* not periodic - extrapolate */

    }

  }

  return;
}
