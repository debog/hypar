/*! @file Euler1DEigen.c
    @author Debojyoti Ghosh
    @brief Contains the functions to compute left and right eigenvectors for the 1D Euler equations.
*/
#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

/*! Compute the left eigenvections for the 1D Euler equations. This function
    just calls the macro #_Euler1DLeftEigenvectors_ and is not used by any
    functions within the 1D Euler module. However, it's necessary to define it
    and provide it to the the solver object (#HyPar) so that it can then send it
    to interpolation functions for a characteristic-based reconstruction.
*/
int Euler1DLeftEigenvectors(
                            double *a_u, /*!< Conserved solution at a grid point */
                            double *a_L, /*!< Array of size nvar^2 = 3^2 to save the matrix of
                                            left eigenvectors in (row-major format). */
                            void   *a_p, /*!< Object of type #Euler1D with physics-related variables */
                            int     a_dir   /*!< Spatial dimension (not used, since this is a 1D system) */
                           )
{
  Euler1D *param  = (Euler1D*)  a_p;
  _Euler1DLeftEigenvectors_(a_u,a_L,param,a_dir);
  return(0);
}

/*! Compute the right eigenvections for the 1D Euler equations. This function
    just calls the macro #_Euler1DRightEigenvectors_ and is not used by any
    functions within the 1D Euler module. However, it's necessary to define it
    and provide it to the the solver object (#HyPar) so that it can then send it
    to interpolation functions for a characteristic-based reconstruction.
*/
int Euler1DRightEigenvectors(
                              double  *a_u, /*!< Conserved solution at a grid point */
                              double  *a_R, /*!< Array of size nvar^2 = 3^2 to save the matrix of
                                               right eigenvectors in (row-major format). */
                              void    *a_p, /*!< Object of type #Euler1D with physics-related variables */
                              int     a_dir   /*!< Spatial dimension (not used, since this is a 1D system) */
                            )
{
  Euler1D *param  = (Euler1D*)  a_p;
  _Euler1DRightEigenvectors_(a_u,a_R,param,a_dir);
  return(0);
}
