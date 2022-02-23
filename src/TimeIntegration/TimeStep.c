/*! @file TimeStep.c
    @brief Advance one time step
    @author Debojyoti Ghosh
*/

#include <basic.h>
#include <simulation_object.h>
#include <timeintegration.h>

/*!
  Advance one time step.
*/
int TimeStep(void *ts /*!< Object of type #TimeIntegration */)
{
  TimeIntegration *TS  = (TimeIntegration*) ts;
  SimulationObject* sim = (SimulationObject*) TS->simulation;

  if (TS->TimeIntegrate) { TS->TimeIntegrate(TS); }

  return(0);
}
