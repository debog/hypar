/*! @file SparseGridsInitialize.cpp
    @brief Initialize all the stuff needed for a sparse grids simulation
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <basic.h>
#include <arrayfunctions.h>
#include <sparse_grids_simulation.h>

/*! Initialize all the stuff required for a sparse grids simulation.
 * Before this function is called, the solver inputs for the full
 * grid object has already been read, and therefore, the global dimensions
 * and processor decompositions are available. This function does the following:
*/
int SparseGridsSimulation::Initialize()
{
  int ierr;

  /* find out the number of spatial dimensions */
  m_ndims = m_sim_fg->solver.ndims;

  /* make sure full grid simulation object is allocated */
  if (m_sim_fg == NULL) {
    fprintf(stderr, "Error in SparseGridsSimulation::Initialize()\n");
    fprintf(stderr, "  m_sim_fg is NULL on rank %d!\n", m_rank);
    return 1;
  }

  /* set full grid processor distribution based on grid size and
   * number of MPI ranks */
  GridDimensions dim_fg(m_ndims);
  for (int d=0; d<m_ndims; d++) dim_fg[d] = m_sim_fg->solver.dim_global[d];
  ProcDistribution iproc_fg;
  ComputeProcessorDistribution(iproc_fg, dim_fg);
  for (int d=0; d<m_ndims; d++) m_sim_fg->mpi.iproc[d] = iproc_fg[d];
  MPIBroadcast_integer( m_sim_fg->mpi.iproc,
                        m_ndims,
                        0,
                        &(m_sim_fg->mpi.world)); CHECKERR(ierr);

  ::WriteInputs( (void*) m_sim_fg, 1, m_rank);
  if (!m_rank) {
    printf("Processor distribution for full grid object: ");
    for (int d=0; d<m_ndims; d++) printf(" %3d", iproc_fg[d]);
    printf("\n");
  }

  if (!m_rank) printf("Initializing sparse grids...\n");
  if (!m_rank) printf("  Number of spatial dimensions: %d\n", m_ndims);

  /* some sanity checks */
  ierr = SanityChecks();
  if (ierr) return ierr;

  /* compute the sizes of grids that go into the combination technique */
  if (!m_rank) {
    printf("  Computing sparse grids dimensions...\n");
  }
  ierr = ComputeSGDimsAndCoeffs();
  if (ierr) return ierr;
  m_nsims_sg = m_combination.size();
  if (!m_rank) {
    printf( "  Number of sparse grid domains in combination technique: %d\n",
            m_nsims_sg );
  }
  if (m_nsims_sg == 0) {
    fprintf(stderr, "Error in SparseGridsSimulation::Initialize()\n");
    fprintf(stderr, "  ComputeSGDimsAndCoeffs() returned empty vector!\n");
    return 1;
  }

  /* compute processor distributions for the sparse grids */
  if (!m_rank) {
    printf("  Computing processor decompositions...\n");
  }
  m_iprocs.resize(m_nsims_sg);
  for (int i = 0; i < m_nsims_sg; i++) {
    ierr = ComputeProcessorDistribution(  m_iprocs[i],
                                          m_combination[i]._dim_ );

    if (ierr) return ierr;
  }

  /* Print some stuff to screen */
  if (!m_rank) {
    printf("Sparse Grids: Combination technique grids sizes and coefficients are:-\n");
    for (int i = 0; i < m_nsims_sg; i++) {
      printf("  %3d: dim = (", i);
      for (int d=0; d<m_ndims; d++) printf(" %4d ", m_combination[i]._dim_[d]);
      printf("),  coeff = %+1.2e, ", m_combination[i]._coeff_);
      printf("iproc = (");
      for (int d=0; d<m_ndims; d++) printf(" %4d ", m_iprocs[i][d]);
      printf(")\n");
    }
  }

  /* allocate full grid data structures (only what is needed) */
  ierr = InitializeBarebones(m_sim_fg);
  if (ierr) return ierr;

  /* create sparse grids simulation objects */
  m_sims_sg.resize(m_nsims_sg);
  /* set the solver parameters for the sparse grids sim objects */
  for (int i = 0; i < m_nsims_sg; i++) {
#ifndef serial
    MPI_Comm_dup(MPI_COMM_WORLD, &(m_sims_sg[i].mpi.world));
#endif
    ierr = SetSolverParameters( m_sims_sg[i],
                                m_combination[i]._dim_,
                                m_iprocs[i],
                                *m_sim_fg,
                                i,
                                m_nsims_sg );
    if (ierr) return ierr;
  }

  /* allocate their data structures (everything) */
  ierr = ::Initialize( (void*) m_sims_sg.data(), m_nsims_sg);
  if (ierr) return ierr;

  /* compute and report total NDOFs for sparse grids */
  long ndof_sg = 0;
  for (int i = 0; i < m_nsims_sg; i++) {
    ndof_sg += m_sims_sg[i].solver.npoints_global;
  }
  long ndof_fg = m_sim_fg->solver.npoints_global;
  if (!m_rank) {
    printf("Total number of DOFs:-\n");
    printf("  using sparse grids: %d\n", (int)ndof_sg);
    printf("  using conventional grid: %d\n", (int)ndof_fg);
  }

  /* done */
  return 0;
}
