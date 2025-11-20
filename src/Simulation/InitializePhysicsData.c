/*! @file InitializePhysicsData.c
    @author Debojyoti Ghosh
    @brief Read in any arrays/data that physics model wants
*/

#include <stdio.h>
#include <simulation_object.h>

/*! For each simulation object, call the physics-specific function to
    read in any physics data that is not a part of the solution vector.
*/
int InitializePhysicsData(void  *a_s,       /*!< Simulation object of type #SimulationObject */
                          int   a_idx,      /*!< Index of this simulation object */
                          int   a_nsims,    /*!< Total number of simuations */
                          int   *a_dim_data /*!< Dimenions of physics-specific data */
                         )
{
  SimulationObject *sim     = (SimulationObject*) a_s;
  HyPar            *solver  = &(sim->solver);
  MPIVariables     *mpi     = &(sim->mpi);

  if (solver->PhysicsInput) {
    int ierr = solver->PhysicsInput(solver, mpi, a_idx, a_nsims, a_dim_data);
    if (ierr) {
      fprintf(stderr, "Error in InitializePhysicsData():\n");
      fprintf(stderr, "  solver->PhysicsInput() returned error %d on rank %d\n",
              ierr, mpi->m_rank);
      return ierr;
    }
  }

  return 0;
}

