/*! @file template_model.h
    @brief [REPLACE] Brief description of your physical model
    @author [YOUR NAME]
    
    [REPLACE] Detailed mathematical description:
    \f{equation}{
      \frac{\partial u}{\partial t} + \nabla \cdot {\bf f}(u) = S(u)
    \f}
    
    References:
    + [REPLACE] Add relevant citations here
*/

#ifndef _TEMPLATE_MODEL_H_
#define _TEMPLATE_MODEL_H_

/*! \def _TEMPLATE_MODEL_
    [REPLACE] Model identifier string - must be unique and match solver.inp
*/
#define _TEMPLATE_MODEL_ "template_model"

#include <basic.h>

/*! \brief Structure containing model-specific parameters and variables
 *
 *  This structure contains the physical parameters, variables, and arrays
 *  specific to your physical model.
*/
typedef struct template_model_parameters {

  /* Physical parameters */
  double param1;              /*!< [REPLACE] Description of parameter 1 */
  double param2;              /*!< [REPLACE] Description of parameter 2 */
  
  /* Model-specific arrays (allocate in Initialize, free in Cleanup) */
  double *field_array;        /*!< [REPLACE] Description of field (e.g., gravity field) */
  int    field_size;          /*!< Size of field_array */
  
  /* Flags and options */
  int    use_feature;         /*!< [REPLACE] Flag to enable/disable feature (0 or 1) */
  char   option_string[_MAX_STRING_SIZE_];  /*!< [REPLACE] String option from physics.inp */
  
  /* Non-dimensionalization factors (if applicable) */
  double reference_length;    /*!< Reference length scale */
  double reference_velocity;  /*!< Reference velocity scale */
  double reference_time;      /*!< Reference time scale */

} TemplateModel;

/* Function declarations - required */
int    TemplateModelInitialize     (void*,void*);
int    TemplateModelCleanup        (void*);

/* Additional function declarations as needed */
/* int    TemplateModelFunctions      (void*,void*,double); */

#endif
