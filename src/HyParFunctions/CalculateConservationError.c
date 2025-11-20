/*! @file CalculateConservationError.c
    @author Debojyoti Ghosh
    @brief Compute the conservation error.
*/
#include <math.h>
#include <basic.h>
#include <mathfunctions.h>
#include <mpivars.h>
#include <hypar.h>

/*! Calculates the error (L2) in conservation by computing the difference between
    the initial volume integral of the solution, and the sum of the current
    volume integral and the time integral of the boundary flux from the start
    of the simulation to the current simulation time.
*/
int CalculateConservationError(
                                void *a_s, /*!< Solver object of type #HyPar */
                                void *a_m  /*!< MPI object of type #MPIVariables */
                              )
{
  HyPar         *solver = (HyPar*) a_s;
  int           v,nvars = solver->m_nvars;
  double        error;

  double base[nvars];
  for (v=0; v<nvars; v++) {
    if (absolute(solver->m_volume_integral_initial[v]) > 1.0)
      base[v] = absolute(solver->m_volume_integral_initial[v]);
    else base[v] = 1.0;
  }

  for (v=0; v<nvars; v++) {
    error =  (solver->m_volume_integral[v]+solver->m_total_boundary_integral[v]-solver->m_volume_integral_initial[v])
           * (solver->m_volume_integral[v]+solver->m_total_boundary_integral[v]-solver->m_volume_integral_initial[v]);
    solver->m_conservation_error[v] = sqrt(error)/base[v];
  }

  return(0);
}
