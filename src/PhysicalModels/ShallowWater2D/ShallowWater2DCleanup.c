/*! @file ShallowWater2DCleanup.c
    @author Debojyoti Ghosh
    @brief Contains the function to clean up the 2D shallow water physics module
*/

#include <stdlib.h>
#include <physicalmodels/shallowwater2d.h>

/*! Function to clean up all physics-related allocations for the 2D shallow water equations */
int ShallowWater2DCleanup(
                   void *a_s /*!< Solver object of type #HyPar */
                  )
{
  ShallowWater2D *param  = (ShallowWater2D*) a_s;
  free(param->m_b);
  return(0);
}
