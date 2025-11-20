/*! @file ShallowWater1DEigen.c
    @author Debojyoti Ghosh
    @brief Contains the functions to compute left and right eigenvectors for the 1D shallow water equations.
*/
#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/shallowwater1d.h>
#include <hypar.h>

/*! Compute the left eigenvections for the 1D shallow water equations. This function
    just calls the macro #_ShallowWater1DLeftEigenvectors_ and is not used by any
    functions within the 1D shallow water module. However, it's necessary to define it
    and provide it to the the solver object (#HyPar) so that it can then send it
    to interpolation functions for a characteristic-based reconstruction.
*/
int ShallowWater1DLeftEigenvectors(
                            double *a_u, /*!< Conserved solution at a grid point */
                            double *a_L, /*!< Array of size nvar^2 = 2^2 to save the matrix of
                                            left eigenvectors in (row-major format). */
                            void   *a_p, /*!< Object of type #ShallowWater1D with physics-related variables */
                            int     a_dir   /*!< Spatial dimension (not used, since this is a 1D system) */
                           )
{
  ShallowWater1D *param  = (ShallowWater1D*)  a_p;
  _ShallowWater1DLeftEigenvectors_(a_u,a_L,param,a_dir);
  return(0);
}

/*! Compute the right eigenvections for the 1D shallow water equations. This function
    just calls the macro #_ShallowWater1DRightEigenvectors_ and is not used by any
    functions within the 1D shallow water module. However, it's necessary to define it
    and provide it to the the solver object (#HyPar) so that it can then send it
    to interpolation functions for a characteristic-based reconstruction.
*/
int ShallowWater1DRightEigenvectors(
                              double  *a_u, /*!< Conserved solution at a grid point */
                              double  *a_R, /*!< Array of size nvar^2 = 2^2 to save the matrix of
                                               right eigenvectors in (row-major format). */
                              void    *a_p, /*!< Object of type #ShallowWater1D with physics-related variables */
                              int     a_dir   /*!< Spatial dimension (not used, since this is a 1D system) */
                            )
{
  ShallowWater1D *param  = (ShallowWater1D*)  a_p;
  _ShallowWater1DRightEigenvectors_(a_u,a_R,param,a_dir);
  return(0);
}
