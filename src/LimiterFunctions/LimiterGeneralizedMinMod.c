/*! @file LimiterGeneralizedMinMod.c
    @author Debojyoti Ghosh
    @brief Generalized MinMod limiter
*/

#include <mathfunctions.h>
#include <limiters.h>

/*! Generalized MinMod limiter

    Reference:
    + Van Leer, B. (1979), "Towards the ultimate conservative difference scheme V.
      A second order sequel to Godunov's method", J. Comput. Phys., 32: 101–136,
      doi:10.1016/0021-9991(79)90145-1
*/
double LimiterGeneralizedMinMod (
                                  double r /*!< Input slope ratio */
                                )
{
  double theta = 1.0;
  double retval = max(0.0,min3(theta*r,0.5*(1.0+r),theta));
  return retval;
}
