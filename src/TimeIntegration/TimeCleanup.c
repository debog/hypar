/*! @file TimeCleanup.c
    @brief Clean up time integration
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(HAVE_CUDA)
#include <arrayfunctions_gpu.h>
#endif
#include <simulation_object.h>
#include <timeintegration.h>

/*!
  Clean up all allocations related to time integration
*/
int TimeCleanup(void *ts /*!< Object of type #TimeIntegration*/)
{
  TimeIntegration* TS = (TimeIntegration*) ts;
  SimulationObject* sim = (SimulationObject*) TS->simulation;
  int ns, nsims = TS->nsims;

  /* close files opened for writing */
  if (!TS->rank) if (sim[0].solver.write_residual) fclose((FILE*)TS->ResidualFile);

#if defined(HAVE_CUDA)
  if (sim[0].solver.use_gpu) {
    if (!strcmp(sim[0].solver.time_scheme,_RK_)) {
      gpuFree(TS->gpu_U);
      gpuFree(TS->gpu_Udot);
      gpuFree(TS->gpu_BoundaryFlux);
    }
  } else {
#endif
    if (!strcmp(sim[0].solver.time_scheme,_RK_)) {
      int i;
      ExplicitRKParameters  *params = (ExplicitRKParameters*)  sim[0].solver.msti;
      for (i=0; i<params->nstages; i++) free(TS->U[i]);            free(TS->U);
      for (i=0; i<params->nstages; i++) free(TS->Udot[i]);         free(TS->Udot);
      for (i=0; i<params->nstages; i++) free(TS->BoundaryFlux[i]); free(TS->BoundaryFlux);
    } else if (!strcmp(sim[0].solver.time_scheme,_FORWARD_EULER_)) {
      int nstages = 1, i;
      for (i=0; i<nstages; i++) free(TS->BoundaryFlux[i]); free(TS->BoundaryFlux);
    } else if (!strcmp(sim[0].solver.time_scheme,_GLM_GEE_)) {
      int i;
      GLMGEEParameters  *params = (GLMGEEParameters*)  sim[0].solver.msti;
      for (i=0; i<2*params->r-1  ; i++) free(TS->U[i]);            free(TS->U);
      for (i=0; i<params->nstages; i++) free(TS->Udot[i]);         free(TS->Udot);
      for (i=0; i<params->nstages; i++) free(TS->BoundaryFlux[i]); free(TS->BoundaryFlux);
    }
#if defined(HAVE_CUDA)
  }
#endif

  /* deallocate arrays */
  free(TS->u_offsets);
  free(TS->u_sizes);
  free(TS->bf_offsets);
  free(TS->bf_sizes);
  free(TS->u  );
  free(TS->rhs);
  for (ns = 0; ns < nsims; ns++) {
    sim[ns].solver.time_integrator = NULL;
  }
  return(0);
}
