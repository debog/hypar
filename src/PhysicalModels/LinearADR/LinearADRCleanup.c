/*! @file LinearADRCleanup.c
    @author Debojyoti Ghosh
    @brief Clean up the linear advection-diffusion-reaction model
*/

#include <stdlib.h>
#include <physicalmodels/linearadr.h>

/*! Clean up the linear advection-diffusion-reaction model */
int LinearADRCleanup(void *s /*!< Solver object of type #HyPar */)
{
  LinearADR *physics = (LinearADR*) s;
  if (physics->a) free(physics->a);
  if (physics->d) free(physics->d);

  return(0);
}
