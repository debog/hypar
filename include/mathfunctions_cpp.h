/*! @file mathfunctions_cpp.h
    @brief Contains function definitions for common mathematical functions for C++ code.
    @author Debojyoti Ghosh
 */

#ifndef _MATHFUNC_CPP_H_
#define _MATHFUNC_CPP_H_

/*! Function to calculate the grid points corresponding to a given interval */
extern "C" void FindInterval(double,double,double*,int,int*,int*);

/*! Fill the ghost cells of a global n-dimensional array */
extern "C" void fillGhostCells( const int* const,
                                const int,
                                double* const,
                                const int,
                                const int,
                                const int* const);

/*! Function to compute trilinear interpolation coefficients */
extern "C" void TrilinearInterpCoeffs(  double,
                                        double,
                                        double,
                                        double,
                                        double,
                                        double,
                                        double,
                                        double,
                                        double,
                                        double*);

/*! Interpolate an n-dimensional grid variable from one grid to another */
extern "C" int InterpolateGlobalnDVar(const int* const,
                                      double** const,
                                      const int* const,
                                      double* const,
                                      const int,
                                      const int,
                                      const int,
                                      const int* const);

#endif
