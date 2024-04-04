/*! @file linearadr.h
    @brief Linear Advection-Diffusion-Reaction model
    @author Debojyoti Ghosh

    Linear Advection-Diffusion-Reachtion\n

    \f{equation}{
      \frac {\partial u} {\partial t}
    + \sum_d \frac {\partial} {\partial x_d} \left( a_d d \right)
    = \sum_d \nu_d \frac {\partial^2 u} {\partial x_d^2}
    + k u
    \f}
    where \f$a_d\f$ are the advection speeds,
    \f$\nu_d\f$ are the diffusion coefficients, and
    \f$k\f$ is the reaction rate.\n
*/

/*! \def _LINEAR_ADVECTION_DIFFUSION_REACTION_
    Linear advection-diffusion-reaction model
*/
#define _LINEAR_ADVECTION_DIFFUSION_REACTION_  "linear-advection-diffusion-reaction"

#include <basic.h>

/*! \def LinearADR
    \brief Structure containing variables and parameters specific
    to the linear advection-diffusion-reaction model.
 *  This structure contains the physical parameters, variables, and function pointers
 *  specific to the linear advection-diffusion-reaction equation.
*/
/*! \brief Structure containing variables and parameters specific
 * to the linear advection-diffusion-reaction model.
 *
 *  This structure contains the physical parameters, variables, and function pointers
 *  specific to the linear advection-diffusion-reaction equation.
*/
typedef struct linearadr_parameters {

  /*! Is the advection field constant (1) or spatially-varying (0) */
  int constant_advection;

  /*! Filename of file to read in spatially-varying advection field from */
  char adv_filename[_MAX_STRING_SIZE_];

  /*! Size of the advection array: depends on whether advection coeff
      is constant or spatially-varying */
  int adv_arr_size;

  /*! advection speed for each variable along each dimension */
  double *a;

  /*! diffusion coefficient for each variable along each dimension */
  double *d;

  /*! turn off upwinding and just take arithmetic average of left
      and right biased fluxes? */
  char centered_flux[_MAX_STRING_SIZE_];

} LinearADR;

/*! Initialize the linear advection-diffusion-reaction model */
int LinearADRInitialize        (void*,void*);
/*! Clean up the linear advection-diffusion-reaction model */
int LinearADRCleanup           (void*);

