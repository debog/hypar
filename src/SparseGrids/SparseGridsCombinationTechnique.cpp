/*! @file SparseGridsCombinationTechnique.cpp
    @brief Functions to implement the combination technique
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <arrayfunctions.h>
#include <std_vec_ops.h>
#include <sparse_grids_simulation.h>

extern "C" void CombineSolutions( SimulationObject*,
                                  double* const *const,
                                  const int,
                                  SimulationObject*,
                                  double* const,
                                  const double* const );

/*! Implements the combination technique where the solutions from all
 *  the sparse grids are combined to give a higher-resolution solution.
 *
 *  The sparse grids domains may have different processor layouts, so this
 *  combination is carried out on rank 0, and the solution is distributed.
*/
void SparseGridsSimulation::CombinationTechnique(SimulationObject* const a_sim /*!< target simulation object on which to combine */)
{
  double** u_sg = (double**) calloc (m_nsims_sg, sizeof(double*));
  std::vector<double> coeffs(m_nsims_sg, 0.0);
  for (int n=0; n<m_nsims_sg; n++) {
    u_sg[n] = m_sims_sg[n].solver.u;
    coeffs[n] = m_combination[n]._coeff_;
  }

  ::CombineSolutions( m_sims_sg.data(),
                      u_sg,
                      m_nsims_sg,
                      a_sim,
                      a_sim->solver.u,
                      coeffs.data() );

  /* done */
  return;
}

