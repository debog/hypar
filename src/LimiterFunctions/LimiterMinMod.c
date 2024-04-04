/*! @file LimiterMinMod.c
    @author Debojyoti Ghosh
    @brief MinMod limiter
*/

#include <mathfunctions.h>
#include <limiters.h>

/*! MinMod limiter

    Reference:
    + Roe, P.L. (1986), "Characteristic-based schemes for the Euler equations",
      Annu. Rev. Fluid Mech., 18: 337–365,
      doi:10.1146/annurev.fl.18.010186.002005
*/
double LimiterMinMod(
                      double r /*!< Input slope ratio */
                    )
{
  double retval = max(0.0,min(1.0,r));
  return retval;
}
