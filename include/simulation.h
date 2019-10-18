/*! @file simulation.h
    @brief Simulation object
    @author Debojyoti Ghosh
*/

#include <mpivars.h>
#include <hypar.h>

/*! \def SimulationObject
 *  \brief Structure defining a simulation
 *  This structure contains an object of type #HyPar
 *  and an object of type #MPIVariables.
*/

/*! \brief Structure defining a simulation
 * 
 * This structure contains an object of type #HyPar
 * and an object of type #MPIVariables.
*/
typedef struct simulation_object {

  MPIVariables  mpi;          /*!< MPI-related variables */
  HyPar         solver;       /*!< Solver-related variables */

} SimulationObject;

/*! Write errors for each simulation */
void SimWriteErrors(void*, int, int, double, double);

int ReadInputs(void*,int,int);/*!< Read the input parameters */
int Initialize(void*,int);/*!< Initialize the solver */
int InitialSolution(void*,int);/*!< Read the initial solution */
int InitializeBoundaries(void*,int);/*!< Initialize the boundary conditions */
int InitializeImmersedBoundaries(void*,int);/*!< Initialize the immersed boundary conditions */
int InitializePhysics(void*,int);/*!< Initialize the physics */
int InitializeSolvers(void*,int);/*!< Initialize the solvers */

int OutputSolution (void*,int);/*!< Write solutions to file */

int Solve(void*,int, int, int);/*!< Solve the PDE - time-integration */
#ifdef with_petsc
int SolvePETSc(void*,int, int, int);  /*!< Solve the PDE using PETSc TS */
#endif

int Cleanup(void*,int);/*!< Clean up: deallocate all arrays and objects */
