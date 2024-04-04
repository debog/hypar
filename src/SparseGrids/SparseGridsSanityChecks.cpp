/*! @file SparseGridsSanityChecks.cpp
    @brief Sanity checks for sparse grids simulations
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <sparse_grids_simulation.h>

/*! Do some basic sanity checks to make sure that a sparse
 * grids simulation can be carried out for this setup.
 *
 * For the full grid:
 *
 * + Check if grid sizes are same along all dimensions.
 * + Check if grid size is a power of 2.
 * + Check if number of MPI ranks are same along all dimensions.
 * + Check if number of MPI ranks is a power of 2.
*/
int SparseGridsSimulation::SanityChecks()
{
  int *dim_global_fg = m_sim_fg->solver.dim_global;
  int *iproc_fg = m_sim_fg->mpi.iproc;

  /* check if grid sizes are same along all dimensions */
  {
    bool flag = true;
    for (int d=0; d<m_ndims; d++) {
      flag = flag && (dim_global_fg[d] == dim_global_fg[0]);
    }
    if (!flag) {
      fprintf(stderr, "Error in SparseGridsSimulation::SanityChecks()\n");
      fprintf(stderr, "  full grid dimensions are not equal in all dimensions.\n");
      return 1;
    }
  }

  /* check if grid size is a power of 2 */
  {
    bool flag = isPowerOfTwo(dim_global_fg[0]);
    if (!flag) {
      fprintf(stderr, "Error in SparseGridsSimulation::SanityChecks()\n");
      fprintf(stderr, "  full grid dimensions are not a power of 2.\n");
      return 1;
    }
  }

  return 0;
}
