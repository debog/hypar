/*! @file NavierStokes2DEigen.c
    @author Debojyoti Ghosh
    @brief Contains the functions to compute left and right eigenvectors for the 2D Navier Stokes equations.
*/
#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes2d.h>
#include <hypar.h>

/*! Compute the left eigenvections for the 2D Navier Stokes equations. This function
    just calls the macro #_NavierStokes2DLeftEigenvectors_ and is not used by any
    functions within the 2D Navier Stokes module. However, it's necessary to define it
    and provide it to the the solver object (#HyPar) so that it can then send it
    to interpolation functions for a characteristic-based reconstruction.
*/
int NavierStokes2DLeftEigenvectors(
                                    double *a_u, /*!< Conserved solution at a grid point */
                                    double *a_L, /*!< Array of size nvar^2 = 4^2 to save the matrix of
                                                    left eigenvectors in (row-major format). */
                                    void   *a_p, /*!< Object of type #NavierStokes2D with physics-related variables */
                                    int    a_dir /*!< Spatial dimension (x or y) */
                                  )
{
  NavierStokes2D *param  = (NavierStokes2D*)  a_p;
  _NavierStokes2DLeftEigenvectors_(a_u,a_L,param->m_gamma,a_dir);
  return(0);
}

/*! Compute the right eigenvections for the 2D Navier Stokes equations. This function
    just calls the macro #_NavierStokes2DRightEigenvectors_ and is not used by any
    functions within the 2D Navier Stokes module. However, it's necessary to define it
    and provide it to the the solver object (#HyPar) so that it can then send it
    to interpolation functions for a characteristic-based reconstruction.
*/
int NavierStokes2DRightEigenvectors(
                                      double  *a_u, /*!< Conserved solution at a grid point */
                                      double  *a_R, /*!< Array of size nvar^2 = 4^2 to save the matrix of
                                                       right eigenvectors in (row-major format). */
                                      void    *a_p, /*!< Object of type #NavierStokes2D with physics-related variables */
                                      int     a_dir /*!< Spatial dimension (x or y) */
                                    )
{
  NavierStokes2D *param  = (NavierStokes2D*)  a_p;
  _NavierStokes2DRightEigenvectors_(a_u,a_R,param->m_gamma,a_dir);
  return(0);
}
