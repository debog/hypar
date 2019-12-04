/*! @file SparseGridsCombinationTechnique.cpp
    @brief Functions to implement the combination technique
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <arrayfunctions.h>
#include <std_vec_ops.h>
#include <sparse_grids_simulation.h>

/*! Implements the combination technique where the solutions from all
 *  the sparse grids are combined to give a higher-resolution solution.
 *
 *  The sparse grids domains may have different processor layouts, so this
 *  combination is carried out on rank 0, and the solution is distributed.
*/
void SparseGridsSimulation::CombinationTechnique(SimulationObject* const a_sim /*!< target simulation object on which to combine */)
{
  /* get the target grid dimensions */
  GridDimensions dim_target;
  StdVecOps::copyFrom(dim_target, a_sim->solver.dim_global, m_ndims);

  /* get number of vector components */
  int nvars = a_sim->solver.nvars;

  /* array size */
  long size = nvars * StdVecOps::product(dim_target);

  /* allocate data array for combined solution */
  double *ug_combined = NULL;
  if (!m_rank) allocateDataArrays(dim_target, nvars, &ug_combined);

  /* for each sparse grids, interpolate onto the target grid dimension
   * and add to the solution */
  for (int n = 0; n < m_nsims_sg; n++) {

    double *u_sg_interpolated = NULL;
    interpolate(dim_target, &u_sg_interpolated, &m_sims_sg[n]);

    if (!m_rank) {
      double coeff = m_combination[n]._coeff_;
      _ArrayAXPY_(u_sg_interpolated, coeff, ug_combined, size);
    }
  }

  /* now partition the global combined solution to its processor */
  MPIPartitionArraynD(  m_ndims,
                        (void*) &(a_sim->mpi),
                        (m_rank ? NULL : ug_combined),
                        a_sim->solver.u,
                        a_sim->solver.dim_global,
                        a_sim->solver.dim_local,
                        a_sim->solver.ghosts,
                        nvars );

  /* free memory */
  if (!m_rank) delete[] ug_combined;

  /* done */
  return;
}

