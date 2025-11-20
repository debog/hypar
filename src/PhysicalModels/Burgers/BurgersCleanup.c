/*! @file BurgersCleanup.c
    @author John Loffeld
    @brief Contains the function to clean up the Burgers physics module
*/

#include <stdlib.h>
#include <physicalmodels/burgers.h>

/*! Function to clean up all physics-related allocations for the Burgers equations */
int BurgersCleanup(void *a_s /*!< Solver object of type #HyPar */)
{
  Burgers *physics = (Burgers*) a_s;

  return(0);
}
