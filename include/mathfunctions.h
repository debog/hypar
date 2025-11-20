/*! @file mathfunctions.h
    @brief Contains function definitions for common mathematical functions.
    @author Debojyoti Ghosh
 */

#ifndef _MATHFUNC_H_
#define _MATHFUNC_H_

#include <math_ops.h>

/*! Function to calculate the grid points corresponding to a given interval */
void FindInterval(double a_val, double a_tolerance, double* a_grid, int a_n, int* a_index_low, int* a_index_high);

/*! Fill the ghost cells of a global n-dimensional array */
void FillGhostCells(const int* const a_dim, const int a_ngpt, double* const a_u, const int a_nvars, const int a_ndims, const int* const a_periodic);

/*! Function to compute trilinear interpolation coefficients */
void TrilinearInterpCoeffs(double a_x, double a_y, double a_z, double a_x1, double a_y1, double a_z1, double a_x2, double a_y2, double a_z2, double* a_coeffs);

/*! Interpolate an n-dimensional grid variable from one grid to another */
int InterpolateGlobalnDVar(const int* const a_dim_src,
                           double** const a_x_src,
                           const int* const a_dim_dst,
                           double* const a_u,
                           const int a_nvars,
                           const int a_ndims,
                           const int a_ngpt,
                           const int* const a_periodic);

#endif
