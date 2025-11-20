/*! @file SparseGridsInitializationWrapup.cpp
    @brief Wrap up all initializations
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <arrayfunctions.h>
#include <sparse_grids_simulation.h>

extern "C" int VolumeIntegral (double*,double*,void*,void*);

/*! Finish all the initializations: Interpolate the initial
    solution on the full grid to all the sparse grids, now
    that the boundary conditions are available.
*/
int SparseGridsSimulation::InitializationWrapup()
{
  int ierr;

  for (int n = 0; n < m_nsims_sg; n++) {

    if (!m_rank) {
      printf("Interpolating initial solution to sparse grids domain %d.\n", n);
    }
    interpolate( &m_sims_sg[n], m_sim_fg );

    HyPar* solver = &(m_sims_sg[n].solver);
    MPIVariables* mpi = &(m_sims_sg[n].mpi);

    int nvars = solver->m_nvars;
    int ghosts = solver->m_ghosts;
    int *dim_local = solver->m_dim_local;

    /* exchange MPI-boundary values of u between processors */
    MPIExchangeBoundariesnD(  m_ndims,
                              nvars,
                              dim_local,
                              ghosts,
                              mpi,
                              solver->m_u  );

    /* calculate volume integral of the initial solution */
    ierr = ::VolumeIntegral(solver->m_volume_integral_initial,
                            solver->m_u,
                            solver,
                            mpi );
    if (ierr) return ierr;

    if (!mpi->m_rank) {
      printf("  Volume integral of the initial solution on sparse grids domain %d:\n", n);
      for (int d=0; d<nvars; d++) {
        printf("    %2d:  %1.16E\n",d,solver->m_volume_integral_initial[d]);
      }
    }

    /* Set initial total boundary flux integral to zero */
    _ArraySetValue_( solver->m_total_boundary_integral, nvars, 0 );

  }

  return 0;
}
