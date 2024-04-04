/*! @file SparseGridsInitialSolution.cpp
    @brief Read in initial solution and intepolate it to the sparse grids.
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <arrayfunctions.h>
#include <sparse_grids_simulation.h>

extern "C" int VolumeIntegral (double*,double*,void*,void*);

/*! Reads in the initial solution for the full grid, saves it in the full
    grid object.  Then interpolates (coarsens) the grid coordinates to all
    the sparge grids. The solution is not yet interpolated because boundary
    information is not yet available. It is done in
    SparseGridsSimulation::InitializationWrapup().
*/
int SparseGridsSimulation::InitialSolution()
{
  int ierr;

  /* first, read in the full grid initial solution */
  int retval = ::InitialSolution( (void*) m_sim_fg, 1 );
  if (retval) return retval;

  /* now, interpolate to all the sparse grids */
  for (int n = 0; n < m_nsims_sg; n++) {

    if (!m_rank) {
      printf("Interpolating grid coordinates to sparse grids domain %d.\n", n);
    }
    interpolateGrid( &m_sims_sg[n], m_sim_fg );

    HyPar* solver = &(m_sims_sg[n].solver);
    MPIVariables* mpi = &(m_sims_sg[n].mpi);

    int nvars = solver->nvars;
    int ghosts = solver->ghosts;
    int *dim_local = solver->dim_local;

    /* calculate dxinv */
    {
      int offset = 0;
      for (int d = 0; d < m_ndims; d++) {
        for (int i = 0; i < dim_local[d]; i++) {
          solver->dxinv[i+offset+ghosts]
            = 2.0 / (   solver->x[i+1+offset+ghosts]
                      - solver->x[i-1+offset+ghosts]  );
        }
        offset += (dim_local[d] + 2*ghosts);
      }
    }

    /* exchange MPI-boundary values of dxinv between processors */
    {
      int offset = 0;
      for (int d = 0; d < m_ndims; d++) {
        ierr = MPIExchangeBoundaries1D( mpi,
                                        &(solver->dxinv[offset]),
                                        dim_local[d],
                                        ghosts,
                                        d,
                                        m_ndims );
        if (ierr) return ierr;
        offset += (dim_local[d] + 2*ghosts);
      }
    }

    /* fill in ghost values of dxinv at physical boundaries by extrapolation */
    {
      int offset = 0;
      for (int d = 0; d < m_ndims; d++) {
        double *dxinv = &(solver->dxinv[offset]);
        int    *dim = dim_local;
        if (mpi->ip[d] == 0) {
          /* fill left boundary along this dimension */
          for (int i = 0; i < ghosts; i++) dxinv[i] = dxinv[ghosts];
        }
        if (mpi->ip[d] == mpi->iproc[d]-1) {
          /* fill right boundary along this dimension */
          for (int i = dim[d]+ghosts; i < dim[d]+2*ghosts; i++) {
            dxinv[i] = dxinv[dim[d]+ghosts-1];
          }
        }
        offset  += (dim[d] + 2*ghosts);
      }
    }

  }

  return 0;
}
