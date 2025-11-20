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
int TimeCleanup(void *a_ts /*!< Object of type #TimeIntegration*/)
{
  TimeIntegration* TS = (TimeIntegration*) a_ts;
  SimulationObject* sim = (SimulationObject*) TS->m_simulation;
  int ns, nsims = TS->m_nsims;

  /* close files opened for writing */
  if (!TS->m_rank) if (sim[0].solver.m_write_residual) fclose((FILE*)TS->m_ResidualFile);

#if defined(HAVE_CUDA)
  if (sim[0].solver.m_use_gpu) {
    if (!strcmp(sim[0].solver.m_time_scheme,_RK_)) {
      gpuFree(TS->m_gpu_U);
      gpuFree(TS->m_gpu_Udot);
      gpuFree(TS->m_gpu_BoundaryFlux);
    }
  } else {
#endif
    if (!strcmp(sim[0].solver.m_time_scheme,_RK_)) {
      int i;
      ExplicitRKParameters  *params = (ExplicitRKParameters*)  sim[0].solver.m_msti;
      for (i=0; i<params->nstages; i++) free(TS->m_U[i]);            free(TS->m_U);
      for (i=0; i<params->nstages; i++) free(TS->m_Udot[i]);         free(TS->m_Udot);
      for (i=0; i<params->nstages; i++) free(TS->m_BoundaryFlux[i]); free(TS->m_BoundaryFlux);
    } else if (!strcmp(sim[0].solver.m_time_scheme,_FORWARD_EULER_)) {
      int nstages = 1, i;
      for (i=0; i<nstages; i++) free(TS->m_BoundaryFlux[i]); free(TS->m_BoundaryFlux);
    } else if (!strcmp(sim[0].solver.m_time_scheme,_GLM_GEE_)) {
      int i;
      GLMGEEParameters  *params = (GLMGEEParameters*)  sim[0].solver.m_msti;
      for (i=0; i<2*params->r-1  ; i++) free(TS->m_U[i]);            free(TS->m_U);
      for (i=0; i<params->nstages; i++) free(TS->m_Udot[i]);         free(TS->m_Udot);
      for (i=0; i<params->nstages; i++) free(TS->m_BoundaryFlux[i]); free(TS->m_BoundaryFlux);
    }
#if defined(HAVE_CUDA)
  }
#endif

  /* deallocate arrays */
  free(TS->m_u_offsets);
  free(TS->m_u_sizes);
  free(TS->m_bf_offsets);
  free(TS->m_bf_sizes);
  free(TS->m_u  );
  free(TS->m_rhs);
  for (ns = 0; ns < nsims; ns++) {
    sim[ns].solver.m_time_integrator = NULL;
  }
  return(0);
}
