/*! @file FindInterval.c
    @author Debojyoti Ghosh
    @brief Find grid point indices corresponding to a spatial interval.
*/

#include <mathfunctions.h>

/*! Given an interval \f$\left[a,b\right], a\leq b\f$, find grid indices \a imin
    and \a imax, such that
    \f{align}{
      imin &= \min\ i\ {\rm satisfying}\ x_i \geq a\\
      imax &= \max\ i\ {\rm satisfying}\  x_i \leq b
    \f}
    where \f$\left\{x_i; 0\leq i < N , x_i < x_{i+1} \forall i \right\}\f$
    represents a 1-dimensional grid.
    \n\n
    Note: This function handles 1-dimensional intervals and grids only.
*/
void FindInterval(
                  double  a_a,      /*!< Lower bound of interval */
                  double  a_b,      /*!< Upper bound of interval */
                  double  *a_x,     /*!< Array of spatial coordinates representing a grid */
                  int     a_N,      /*!< Number of grid points / size of x */
                  int     *a_imin,  /*!< Lowest grid index within [a_a,b] */
                  int     *a_imax   /*!< Highest grid index within [a_a,b] */
                 )
{
  int i;
  *a_imax = -1;
  *a_imin =  a_N;

  double min_dx = a_x[1] - a_x[0];
  for (i = 2; i < a_N; i++) {
    double dx = a_x[i] - a_x[i-1];
    if (dx < min_dx) min_dx = dx;
  }
  double tol = 1e-10 * min_dx;

  for (i = 0; i < a_N; i++) {
    if (a_x[i] <= (a_b+tol)) *a_imax = i+1;
  }
  for (i = a_N-1; i > -1; i--) {
    if (a_x[i] >= (a_a-tol)) *a_imin = i;
  }

  return;
}
