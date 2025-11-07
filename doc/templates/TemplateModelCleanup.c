/*! @file TemplateModelCleanup.c
    @author [YOUR NAME]
    @brief Clean up the Template Model physics module
*/

#include <stdlib.h>
#include <physicalmodels/template_model.h>

/*!
 * Function to clean up all physics-related allocations for the Template Model.
 *
 * This function is called at the end of the simulation to free any memory
 * that was allocated in the Initialize function or during the simulation.
*/
int TemplateModelCleanup(void *s /*!< Object of type TemplateModel (cast from void*) */)
{
  TemplateModel *physics = (TemplateModel*) s;

  /* [REPLACE] Free any allocated arrays */

  /* Example: */
  /*
  if (physics->field_array) {
    free(physics->field_array);
    physics->field_array = NULL;
  }
  */

  return(0);
}
