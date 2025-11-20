/*! @file ShallowWater2DFunctions.c
    @author Debojyoti Ghosh
    @brief Miscellaneous functions for the 2D shallow water equations.
*/

#include <math.h>
#include <basic.h>
#include <physicalmodels/shallowwater2d.h>

/*! Compute the Roe-averaged state for the 2D shallow water equations. This function
    just calls the macro #_ShallowWater2DRoeAverage_ and is not used by any
    functions within the 2D shallow water module. However, it's necessary to define it
    and provide it to the the solver object (#HyPar) so that it can then send it
    to interpolation functions for a characteristic-based reconstruction.
*/
int ShallowWater2DRoeAverage(
                      double  *a_uavg, /*!< The computed Roe-averaged state */
                      double  *a_uL,   /*!< Left state (conserved variables)*/
                      double  *a_uR,   /*!< Right state (conserved variables)*/
                      void    *a_p     /*!< Object of type #ShallowWater2D with physics-related variables */
                     )
{
  ShallowWater2D *param  = (ShallowWater2D*) a_p;
  _ShallowWater2DRoeAverage_(a_uavg,a_uL,a_uR,param);
  return(0);
}
