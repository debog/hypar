/*! @file VlasovCleanup.c
    @author John Loffeld
    @brief Contains the function to clean up the Vlasov physics module
*/

#include <stdlib.h>
#include <physicalmodels/vlasov.h>

/*! Function to clean up all physics-related allocations for the Vlasov equations */
int VlasovCleanup(void *s /*!< Solver object of type #HyPar */)
{
  Vlasov *physics = (Vlasov*) s;
  
  free(physics->e_field);

#ifdef fftw
  free(physics->sum_buffer);

  fftw_destroy_plan(physics->plan_forward);
  fftw_destroy_plan(physics->plan_backward);

  fftw_free(physics->phys_buffer);
  fftw_free(physics->fourier_buffer);
#endif

  return(0);
}
