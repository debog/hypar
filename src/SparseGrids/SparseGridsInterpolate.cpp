/*! @file SparseGridsInterpolate.cpp
    @brief Functions to interpolate from one grid to another
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <arrayfunctions.h>
#include <mathfunctions_cpp.h>
#include <std_vec_ops.h>
#include <sparse_grids_simulation.h>

/*! Interpolate data from one grid to another. Note that along each dimension,
    the ratio of the number of grid points in the source grid and that in the
    destination grid must be an integer power of 2 (negative or positive).
*/
void SparseGridsSimulation::interpolate(  SimulationObject* const       a_dst,  /*!< Destination object */
                                          const SimulationObject* const a_src   /*!< Source object */
                                      )
{
  /* get the destination grid dimensions */
  GridDimensions dim_dst;
  StdVecOps::copyFrom(dim_dst, a_dst->solver.dim_global, m_ndims);

  /* get number of vector components */
  int nvars = a_src->solver.nvars;
  if (nvars != a_dst->solver.nvars) {
    fprintf(stderr, "Error in SparseGridsSimulation::interpolate(): unequal nvars\n");
    exit(1);
  }

  double *ug_dst = NULL;
  interpolate(dim_dst, &ug_dst, a_src);

  /* partition the destination */
  MPIPartitionArraynDwGhosts( m_ndims,
                              (void*) &(a_dst->mpi),
                              (a_dst->mpi.rank ? NULL : ug_dst),
                              a_dst->solver.u,
                              a_dst->solver.dim_global,
                              a_dst->solver.dim_local,
                              a_dst->solver.ghosts,
                              nvars);

  if (!m_rank) {
    free(ug_dst);
  }

  return;
}

/*! Interpolate data from one grid to another of a desired resolution.
    Note that along each dimension, the ratio of the number of grid points
    in the source grid and that in the destination grid must be an integer
    power of 2 (negative or positive).

    The incoming pointer must be NULL. After this function is executed, it
    will point to a chunk of memory with the interpolated solution. This is
    the *global* solution.

    Since the source and destination grids may have different processor layouts,
    the interpolation is carried out on rank 0 by allocating global arrays on it.
    Therefore, this function is *not* scalable.
*/
void SparseGridsSimulation::interpolate(const GridDimensions& a_dim_dst, /*!< grid dimensions to interpolate to */
                                        double** const a_u_dst, /*!< pointer to array containing interpolated data */
                                        const SimulationObject* const a_src /*!< Source object */
                                       )
{
  if ((*a_u_dst) != NULL) {
    fprintf(stderr, "Error in SparseGridsSimulation::interpolate() - \n");
    fprintf(stderr, "  a_u_dst is not NULL!\n");
    exit(1);
  }

  /* get the source grid dimensions */
  GridDimensions dim_src;
  StdVecOps::copyFrom(dim_src, a_src->solver.dim_global, m_ndims);

  /* get number of vector components and number of ghost points*/
  int nvars = a_src->solver.nvars;
  int ghosts = a_src->solver.ghosts;

  /* gather the source on rank 0 */
  double *ug_src = NULL;
  if (!m_rank) allocateDataArrays(dim_src, nvars, &ug_src, ghosts);
  MPIGatherArraynDwGhosts( m_ndims,
                           (void*) &(a_src->mpi),
                           ug_src,
                           a_src->solver.u,
                           a_src->solver.dim_global,
                           a_src->solver.dim_local,
                           ghosts,
                           nvars );

  std::vector<int> periodic_arr(m_ndims);
  for (int i=0; i<m_ndims; i++) {
    periodic_arr[i] = (m_is_periodic[i] ? 1 : 0);
  }

  if (!m_rank) {
    int ierr = ::InterpolateGlobalnDVar(  a_dim_dst.data(),
                                          a_u_dst,
                                          dim_src.data(),
                                          ug_src,
                                          nvars,
                                          ghosts,
                                          m_ndims,
                                          periodic_arr.data());
    if (ierr) {
      fprintf(stderr,"InterpolateGlobalnDVar() returned with error!\n");
      exit(1);
    }
  }

  return;
}
