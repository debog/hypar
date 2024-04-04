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
int TimePostStep(void *ts /*!< Object of type #TimeIntegration */)
{
  TimeIntegration* TS = (TimeIntegration*) ts;
  SimulationObject* sim = (SimulationObject*) TS->simulation;
  int ns, nsims = TS->nsims;

  /* update current time */
  TS->waqt += TS->dt;

  if ((TS->iter+1)%sim[0].solver.screen_op_iter == 0) {

#if defined(HAVE_CUDA)
    if (sim[0].solver.use_gpu) {
      TS->norm = -1;
    } else {
#endif
      /* Calculate norm for this time step */
      double sum = 0.0;
      double npts = 0;
      for (ns = 0; ns < nsims; ns++) {
        _ArrayAXPY_(sim[ns].solver.u,-1.0,(TS->u+TS->u_offsets[ns]),TS->u_sizes[ns]);
        sum += ArraySumSquarenD( sim[ns].solver.nvars,
                                 sim[ns].solver.ndims,
                                 sim[ns].solver.dim_local,
                                 sim[ns].solver.ghosts,
                                 sim[ns].solver.index,
                                 (TS->u+TS->u_offsets[ns]) );
        npts += (double)sim[ns].solver.npoints_global;
      }

      double global_sum = 0;
      MPISum_double(  &global_sum,
                      &sum,1,
                      &(sim[0].mpi.world) );

      TS->norm = sqrt(global_sum/npts);
#if defined(HAVE_CUDA)
    }
#endif

    /* write to file */
    if (TS->ResidualFile) {
      fprintf((FILE*)TS->ResidualFile,"%10d\t%E\t%E\n",TS->iter+1,TS->waqt,TS->norm);
    }

  }


#if defined(HAVE_CUDA)
  if (!sim[0].solver.use_gpu) {
#endif
    for (ns = 0; ns < nsims; ns++) {

      if (!strcmp(sim[ns].solver.ConservationCheck,"yes")) {
        /* calculate volume integral of the solution at this time step */
        IERR sim[ns].solver.VolumeIntegralFunction( sim[ns].solver.VolumeIntegral,
                                                    sim[ns].solver.u,
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
        sim[ns].solver.PostStep(sim[ns].solver.u,&(sim[ns].solver),&(sim[ns].mpi),TS->waqt,TS->iter);
      }

    }
#if defined(HAVE_CUDA)
  }
#endif

  gettimeofday(&TS->iter_end_time,NULL);
  long long walltime;
  walltime = (  (TS->iter_end_time.tv_sec * 1000000 + TS->iter_end_time.tv_usec)
              - (TS->iter_start_time.tv_sec * 1000000 + TS->iter_start_time.tv_usec));
  TS->iter_wctime = (double) walltime / 1000000.0;
  TS->iter_wctime_total += TS->iter_wctime;

  double global_total = 0, global_wctime = 0, global_mpi_total = 0, global_mpi_wctime = 0;

  MPIMax_double(&global_wctime, &TS->iter_wctime, 1, &(sim[0].mpi.world));
  MPIMax_double(&global_total, &TS->iter_wctime_total, 1, &(sim[0].mpi.world));

#if defined(HAVE_CUDA)
  MPIMax_double(&global_mpi_wctime, &sim[0].mpi.wctime, 1, &(sim[0].mpi.world));
  MPIMax_double(&global_mpi_total, &sim[0].mpi.wctime_total, 1, &(sim[0].mpi.world));
#endif

  return(0);
}

