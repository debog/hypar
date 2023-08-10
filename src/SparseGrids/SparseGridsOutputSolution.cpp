/*! @file SparseGridsOutputSolution.cpp
    @brief Output the solution to file
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <sparse_grids_simulation.h>

int OutputSolution (void*,int, double);

/*! Write solutions to file */
void SparseGridsSimulation::OutputSolution(double a_time /*!< simulation time */)
{
  /* if asked for, write individual sparse grids solutions */
  if (m_write_sg_solutions == 1) {
    for (int ns = 0; ns < m_nsims_sg; ns++) {
      if (m_sims_sg[ns].solver.PhysicsOutput) {
        m_sims_sg[ns].solver.PhysicsOutput( &(m_sims_sg[ns].solver),
                                            &(m_sims_sg[ns].mpi),
                                            a_time );
      }
    }
    ::OutputSolution((void*)m_sims_sg.data(), m_nsims_sg, a_time);
  }

  /* Combine the sparse grids solutions to full grid */
  CombinationTechnique(m_sim_fg);

  /* Write the full grid solution */
  ::OutputSolution((void*)m_sim_fg, 1, a_time);

  return;
}

