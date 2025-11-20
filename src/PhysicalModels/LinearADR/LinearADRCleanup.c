/*! @file LinearADRCleanup.c
    @author Debojyoti Ghosh
    @brief Clean up the linear advection-diffusion-reaction model
*/

#include <stdlib.h>
#include <physicalmodels/linearadr.h>

/*! Clean up the linear advection-diffusion-reaction model */
int LinearADRCleanup(void *a_s /*!< Solver object of type #HyPar */)
{
  LinearADR *physics = (LinearADR*) a_s;
  if (physics->m_a) free(physics->m_a);
  if (physics->m_d) free(physics->m_d);

  return(0);
}
