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
  free(physics->potential);

#ifdef fftw
  if(physics->self_consistent_electric_field) {
    free(physics->sum_buffer);

    fftw_destroy_plan(physics->plan_forward_e);
    fftw_destroy_plan(physics->plan_backward_e);

    fftw_free(physics->phys_buffer_e);
    fftw_free(physics->fourier_buffer_e);

    fftw_destroy_plan(physics->plan_forward_phi);
    fftw_destroy_plan(physics->plan_backward_phi);

    fftw_free(physics->phys_buffer_phi);
    fftw_free(physics->fourier_buffer_phi);
  }
#endif

  return(0);
}
