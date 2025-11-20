/*! @file LinearADRJacobian.c
    @author Debojyoti Ghosh
    @brief Function to compute the hyperbolic flux Jacobian for the linear advection-diffusion-reaction system.
*/

#include <mathfunctions.h>
#include <physicalmodels/linearadr.h>

/*! Function to compute the flux Jacobian for the hyperbolic (advection) part of the
    linear-advection-diffusion-reaction model. */
int LinearADRAdvectionJacobian( double* a_Jac,  /*!< Jacobian matrix of size 1 (nvar = 1) */
                                double* a_u,    /*!< solution at a grid point */
                                void*   a_p,    /*!< object containing physics-related parameters */
                                int     a_dir,  /*!< dimension (x/y/z) */
                                int     a_nvars,/*!< number of components */
                                int     a_upw   /*!< 0 -> send back complete Jacobian,
                                                   1 -> send back Jacobian of right(+)-moving flux,
                                                  -1 -> send back Jacobian of left(-)-moving flux*/ )
{
  LinearADR *param = (LinearADR*) a_p;

  if (param->m_a) {
    *a_Jac =    (1-absolute(a_upw))*absolute(param->m_a[a_dir])
           +  absolute(a_upw) * (1+a_upw) * max(0,param->m_a[a_dir]) * 0.5
           -  absolute(a_upw) * (1-a_upw) * min(0,param->m_a[a_dir]) * 0.5 ;
  } else {
    /* no advection term */
    *a_Jac = 0.0;
  }

  return 0;
}

/*! Function to compute the  Jacobian for the parabolic (diffusion) part of the
    linear-advection-diffusion-reaction model. */
int LinearADRDiffusionJacobian( double* a_Jac,  /*!< Jacobian matrix of size 1 (nvar = 1) */
                                double* a_u,    /*!< solution at a grid point */
                                void*   a_p,    /*!< object containing physics-related parameters */
                                int     a_dir,  /*!< dimension (x/y/z) */
                                int     a_nvars /*!< number of components */ )
{
  LinearADR *param = (LinearADR*) a_p;

  int v;
  for (v = 0; v < a_nvars; v++) {
    a_Jac[a_nvars*v+v] = -param->m_d[a_nvars*a_dir+v];
  }

  return 0;
}
