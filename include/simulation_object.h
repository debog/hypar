/*! @file simulation_object.h
    @brief Simulation object
    @author Debojyoti Ghosh
*/

#ifndef _simulation_struct_h_
#define _simulation_struct_h_

#include <mpivars_struct.h>
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
  int           is_barebones; /*!< Indicate whether this simulation object
                                   is only a barebones one */

} SimulationObject;

#endif
