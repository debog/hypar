/*! @file math_ops.h
    @brief Contains macros for common mathematical functions.
    @author Debojyoti Ghosh
 */

#ifndef _MATH_OPS_H_
#define _MATH_OPS_H_

#include <math.h>

/*! \def min
 *  Minimum of two numbers
*/
#define min(a,b) ((a)<(b)?(a):(b))
/*! \def max
 * Maximum of two numbers
*/
#define max(a,b) ((a)>(b)?(a):(b))

/*! \def min3
 *  Minimum of three numbers
*/
#define min3(a,b,c) min(min((a),(b)),min((b),(c)))
/*! \def max3
 *  Maximum of three numbers
*/
#define max3(a,b,c) max(max((a),(b)),max((b),(c)))

/*! \def absolute
 * Absolute value
*/
#define absolute(a) ((a)<0?-(a):(a))

/*! \def raiseto
 * Raise to a power: y = x^a
*/
#define raiseto(x,a) (exp((a)*log(x)))

/*! \def raiseto_int
 * Raise to a power (int only): y = x^a
*/
#define raiseto_int(y,x,a) \
  { \
    int arraycounter; \
    y = x; \
    for (arraycounter=1; arraycounter<a; arraycounter++) { \
      y *= x; \
    } \
  }

/*! \def sign
 * Returns the sign of the argument
*/
#define sign(a) ((a)<0?-1.0:1.0)

#endif
