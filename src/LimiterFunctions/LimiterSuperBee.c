/*! @file LimiterSuperbee.c
    @author Debojyoti Ghosh
    @brief Superbee limiter
*/

#include <mathfunctions.h>
#include <limiters.h>

/*! Superbee limiter

    Reference:
    + Roe, P.L. (1986), "Characteristic-based schemes for the Euler equations",
      Annu. Rev. Fluid Mech., 18: 337–365,
      doi:10.1146/annurev.fl.18.010186.002005
*/
double LimiterSuperBee(
                        double r /*!< Input slope ratio */
                      )
{
  double retval = max3(0, min(2*r,1), min(r,2));
  return retval;
}
