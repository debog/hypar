/*! @file TimePostStep.c
    @brief Post-time-step function
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#if defined(HAVE_CUDA)
#include <arrayfunctions_gpu.h>
#else
#include <arrayfunctions.h>
#endif
#include <mpivars.h>
#include <simulation_object.h>
#include <timeintegration.h>

/*!
  Post-time-step function: this function is called at the end of
  each time step.
  + It updates the current simulation time.
  + It calls functions to print information and to write
    transient solution to file.
  + It will also call any physics-specific post-time-step function,
    if defined.
*/
int TimePostStep(void *a_ts /*!< Object of type #TimeIntegration */)
{
  TimeIntegration* TS = (TimeIntegration*) a_ts;
  SimulationObject* sim = (SimulationObject*) TS->m_simulation;
  int ns, nsims = TS->m_nsims;

  /* update current time */
  TS->m_waqt += TS->m_dt;

  if ((TS->m_iter+1)%sim[0].solver.m_screen_op_iter == 0) {

#if defined(HAVE_CUDA)
    if (sim[0].solver.m_use_gpu) {
      TS->m_norm = -1;
    } else {
#endif
      /* Calculate norm for this time step */
      double sum = 0.0;
      double npts = 0;
      for (ns = 0; ns < nsims; ns++) {
        _ArrayAXPY_(sim[ns].solver.m_u,-1.0,(TS->m_u+TS->m_u_offsets[ns]),TS->m_u_sizes[ns]);
        sum += ArraySumSquarenD( sim[ns].solver.m_nvars,
                                 sim[ns].solver.m_ndims,
                                 sim[ns].solver.m_dim_local,
                                 sim[ns].solver.m_ghosts,
                                 sim[ns].solver.m_index,
                                 (TS->m_u+TS->m_u_offsets[ns]) );
        npts += (double)sim[ns].solver.m_npoints_global;
      }

      double global_sum = 0;
      MPISum_double(  &global_sum,
                      &sum,1,
                      &(sim[0].mpi.m_world) );

      TS->m_norm = sqrt(global_sum/npts);
#if defined(HAVE_CUDA)
    }
#endif

    /* write to file */
    if (TS->m_ResidualFile) {
      fprintf((FILE*)TS->m_ResidualFile,"%10d\t%E\t%E\n",TS->m_iter+1,TS->m_waqt,TS->m_norm);
    }

  }


#if defined(HAVE_CUDA)
  if (!sim[0].solver.m_use_gpu) {
#endif
    for (ns = 0; ns < nsims; ns++) {

      if (!strcmp(sim[ns].solver.m_conservation_check,"yes")) {
        /* calculate volume integral of the solution at this time step */
        IERR sim[ns].solver.VolumeIntegralFunction( sim[ns].solver.m_volume_integral,
                                                    sim[ns].solver.m_u,
                                                    &(sim[ns].solver),
                                                    &(sim[ns].mpi) ); CHECKERR(ierr);
        /* calculate surface integral of the flux at this time step */
        IERR sim[ns].solver.BoundaryIntegralFunction( &(sim[ns].solver),
                                                      &(sim[ns].mpi)); CHECKERR(ierr);
        /* calculate the conservation error at this time step       */
        IERR sim[ns].solver.CalculateConservationError( &(sim[ns].solver),
                                                        &(sim[ns].mpi)); CHECKERR(ierr);
      }

      if (sim[ns].solver.PostStep) {
        sim[ns].solver.PostStep(sim[ns].solver.m_u,&(sim[ns].solver),&(sim[ns].mpi),TS->m_waqt,TS->m_iter);
      }

    }
#if defined(HAVE_CUDA)
  }
#endif

  gettimeofday(&TS->m_iter_end_time,NULL);
  long long walltime;
  walltime = (  (TS->m_iter_end_time.tv_sec * 1000000 + TS->m_iter_end_time.tv_usec)
              - (TS->m_iter_start_time.tv_sec * 1000000 + TS->m_iter_start_time.tv_usec));
  TS->m_iter_wctime = (double) walltime / 1000000.0;
  TS->m_iter_wctime_total += TS->m_iter_wctime;

  double global_total = 0, global_wctime = 0, global_mpi_total = 0, global_mpi_wctime = 0;

  MPIMax_double(&global_wctime, &TS->m_iter_wctime, 1, &(sim[0].mpi.m_world));
  MPIMax_double(&global_total, &TS->m_iter_wctime_total, 1, &(sim[0].mpi.m_world));

#if defined(HAVE_CUDA)
  MPIMax_double(&global_mpi_wctime, &sim[0].mpi.m_wctime, 1, &(sim[0].mpi.m_world));
  MPIMax_double(&global_mpi_total, &sim[0].mpi.m_wctime_total, 1, &(sim[0].mpi.m_world));
#endif

  return(0);
}

