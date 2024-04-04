/*! @file LimiterVanLeer.c
    @author Debojyoti Ghosh
    @brief van Leer limiter
*/

#include <mathfunctions.h>
#include <limiters.h>

/*! van Leer limiter

    Reference:
    + Van Leer, B. (1974), "Towards the ultimate conservative difference scheme II.
      Monotonicity and conservation combined in a second order scheme",
      J. Comput. Phys., 14 (4): 361–370,
      doi:10.1016/0021-9991(74)90019-9
*/
double LimiterVanLeer (
                        double r /*!< Input slope ratio */
                      )
{
  double retval = (r+absolute(r)) / (1.0+absolute(r));
  return retval;
}
