/*! @file SingleSimulationDefine.cpp
    @brief Define a single simulation object
    @author Debojyoti Ghosh
*/

#include <string>
#include <single_simulation.h>

/*! Define the single simulation object */
int SingleSimulation::define( int a_rank, /*!< MPI rank of this process */
                              int a_nproc /*!< Total number of MPI ranks */
                            )
{
  if (m_is_defined) {
    fprintf(stderr,"Error: object already defined on rank %d.\n", a_rank);
    return 1;
  }

  m_rank = a_rank;
  m_nproc = a_nproc;

  m_sim = new SimulationObject;
  m_sim->solver.my_idx = 0;
  m_sim->solver.nsims = 1;
  m_sim->mpi.rank = m_rank;
  m_sim->mpi.nproc = m_nproc;

  if (!m_rank) {
    printf("Allocated simulation object(s).\n");
  }

  m_is_defined = true;
  return 0;
}
