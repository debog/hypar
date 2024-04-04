/*! @file LinearADRJacobian.c
    @author Debojyoti Ghosh
    @brief Function to compute the hyperbolic flux Jacobian for the linear advection-diffusion-reaction system.
*/

#include <mathfunctions.h>
#include <physicalmodels/linearadr.h>

/*! Function to compute the flux Jacobian for the hyperbolic (advection) part of the
    linear-advection-diffusion-reaction model. */
int LinearADRAdvectionJacobian( double* Jac,  /*!< Jacobian matrix of size 1 (nvar = 1) */
                                double* u,    /*!< solution at a grid point */
                                void*   p,    /*!< object containing physics-related parameters */
                                int     dir,  /*!< dimension (x/y/z) */
                                int     nvars,/*!< number of components */
                                int     upw   /*!< 0 -> send back complete Jacobian,
                                                   1 -> send back Jacobian of right(+)-moving flux,
                                                  -1 -> send back Jacobian of left(-)-moving flux*/ )
{
  LinearADR *param = (LinearADR*) p;

  if (param->a) {
    *Jac =    (1-absolute(upw))*absolute(param->a[dir])
           +  absolute(upw) * (1+upw) * max(0,param->a[dir]) * 0.5
           -  absolute(upw) * (1-upw) * min(0,param->a[dir]) * 0.5 ;
  } else {
    /* no advection term */
    *Jac = 0.0;
  }

  return 0;
}

/*! Function to compute the  Jacobian for the parabolic (diffusion) part of the
    linear-advection-diffusion-reaction model. */
int LinearADRDiffusionJacobian( double* Jac,  /*!< Jacobian matrix of size 1 (nvar = 1) */
                                double* u,    /*!< solution at a grid point */
                                void*   p,    /*!< object containing physics-related parameters */
                                int     dir,  /*!< dimension (x/y/z) */
                                int     nvars /*!< number of components */ )
{
  LinearADR *param = (LinearADR*) p;

  int v;
  for (v = 0; v < nvars; v++) {
    Jac[nvars*v+v] = -param->d[nvars*dir+v];
  }

  return 0;
}
