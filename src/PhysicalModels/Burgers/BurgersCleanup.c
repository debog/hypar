/*! @file BurgersCleanup.c
    @author John Loffeld
    @brief Contains the function to clean up the Burgers physics module
*/

#include <stdlib.h>
#include <physicalmodels/burgers.h>

/*! Function to clean up all physics-related allocations for the Burgers equations */
int BurgersCleanup(void *s /*!< Solver object of type #HyPar */)
{
  Burgers *physics = (Burgers*) s;

  return(0);
}
