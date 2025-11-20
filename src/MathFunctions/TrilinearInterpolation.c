/*! @file TrilinearInterpolation.c
    @brief Compute coefficients for trilinear interpolation
    @author Debojyoti Ghosh
*/

#include <mathfunctions.h>

/*!
  This function computes the coefficients for a trilinear interpolation at a given
  point (x,y,z) inside a cube defined by [xmin,xmax] X [ymin,ymax] X [zmin,zmax].
  The coefficients are stored in an array of size 8 with each element corresponding
  to a corner of the cube in the following order:\n
  coeffs[0] => xmin,ymin,zmin\n
  coeffs[1] => xmax,ymin,zmin\n
  coeffs[2] => xmin,ymax,zmin\n
  coeffs[3] => xmax,ymax,zmin\n
  coeffs[4] => xmin,ymin,zmax\n
  coeffs[5] => xmax,ymin,zmax\n
  coeffs[6] => xmin,ymax,zmax\n
  coeffs[7] => xmax,ymax,zmax
*/
void TrilinearInterpCoeffs(
                            double a_xmin,  /*!< x-coordinate of the lower-end */
                            double a_xmax,  /*!< x-coordinate of the higher-end */
                            double a_ymin,  /*!< y-coordinate of the lower-end */
                            double a_ymax,  /*!< y-coordinate of the higher-end */
                            double a_zmin,  /*!< z-coordinate of the lower-end */
                            double a_zmax,  /*!< z-coordinate of the higher-end */
                            double a_x,     /*!< x-coordinate of the point to interpolate at */
                            double a_y,     /*!< y-coordinate of the point to interpolate at */
                            double a_z,     /*!< z-coordinate of the point to interpolate at */
                            double *a_coeffs/*!< array of size 8 (pre-allocated) to store the coefficients in */
                          )
{
  double vol_inv = 1 / ((a_xmax-a_xmin)*(a_ymax-a_ymin)*(a_zmax-a_zmin));
  double tldx1 = a_x - a_xmin;
  double tldx2 = a_xmax - a_x;
  double tldy1 = a_y - a_ymin;
  double tldy2 = a_ymax - a_y;
  double tldz1 = a_z - a_zmin;
  double tldz2 = a_zmax - a_z;

  a_coeffs[0] = tldz2 * tldy2 * tldx2 * vol_inv;
  a_coeffs[1] = tldz2 * tldy2 * tldx1 * vol_inv;
  a_coeffs[2] = tldz2 * tldy1 * tldx2 * vol_inv;
  a_coeffs[3] = tldz2 * tldy1 * tldx1 * vol_inv;
  a_coeffs[4] = tldz1 * tldy2 * tldx2 * vol_inv;
  a_coeffs[5] = tldz1 * tldy2 * tldx1 * vol_inv;
  a_coeffs[6] = tldz1 * tldy1 * tldx2 * vol_inv;
  a_coeffs[7] = tldz1 * tldy1 * tldx1 * vol_inv;

  return;
}
