/*! @file TimeRK.c
    @brief Explicit Runge-Kutta method
    @author Debojyoti Ghosh, Youngdae Kim
*/

#include <basic.h>
#if defined (HAVE_CUDA)
#include <arrayfunctions_gpu.h>
#else
#include <arrayfunctions.h>
#endif
#include <simulation_object.h>
#include <timeintegration.h>
#include <time.h>
#include <math.h>

/*!
  Advance the ODE given by
  \f{equation}{
    \frac{d{\bf u}}{dt} = {\bf F} \left({\bf u}\right)
  \f}
  by one time step of size #HyPar::dt using the forward Euler method
  given by
  \f{align}{
    {\bf U}^{\left(i\right)} &= {\bf u}_n + \Delta t \sum_{j=1}^{i-1} a_{ij} {\bf F}\left({\bf U}^{\left(j\right)}\right), \\
    {\bf u}_{n+1} &= {\bf u}_n + \Delta t \sum_{i=1}^s b_{i} {\bf F}\left({\bf U}^{\left(i\right)}\right),
  \f}
  where the subscript represents the time level, the superscripts represent the stages, \f$\Delta t\f$ is the
  time step size #HyPar::dt, and \f${\bf F}\left({\bf u}\right)\f$ is computed by #TimeIntegration::RHSFunction.
  The Butcher tableaux coefficients are \f$a_{ij}\f$ (#ExplicitRKParameters::A) and \f$b_i\f$
  (#ExplicitRKParameters::b).

  Note: In the code #TimeIntegration::Udot is equivalent to \f${\bf F}\left({\bf u}\right)\f$.
*/
int TimeRK(void *ts /*!< Object of type #TimeIntegration */)
{
  TimeIntegration* TS = (TimeIntegration*) ts;
  SimulationObject* sim = (SimulationObject*) TS->simulation;
  ExplicitRKParameters *params = (ExplicitRKParameters*) sim[0].solver.msti;
  int ns, stage, i, nsims = TS->nsims;

#if defined(HAVE_CUDA)
  if (sim[0].solver.use_gpu) {

    /* Calculate stage values */
    for (stage = 0; stage < params->nstages; stage++) {

      double stagetime = TS->waqt + params->c[stage]*TS->dt;

      for (ns = 0; ns < nsims; ns++) {
          gpuArrayCopy1D(  sim[ns].solver.gpu_u,
                           (TS->gpu_U + stage*TS->u_size_total + TS->u_offsets[ns]),
                           (TS->u_sizes[ns])  );

      }

      for (i = 0; i < stage; i++) {
          gpuArrayAXPY(  TS->gpu_Udot + i*TS->u_size_total,
                         (TS->dt * params->A[stage*params->nstages+i]),
                         TS->gpu_U + stage*TS->u_size_total,
                         TS->u_size_total  );

      }

      for (ns = 0; ns < nsims; ns++) {
        if (sim[ns].solver.PreStage) {
          fprintf(stderr,"ERROR in TimeRK(): Call to solver->PreStage() commented out!\n");
          return 1;
//          sim[ns].solver.PreStage( stage,
//                                   (TS->gpu_U),
//                                   &(sim[ns].solver),
//                                   &(sim[ns].mpi),
//                                   stagetime ); CHECKERR(ierr);
        }
      }

      for (ns = 0; ns < nsims; ns++) {
        if (sim[ns].solver.PostStage) {
          sim[ns].solver.PostStage(  (TS->gpu_U + stage*TS->u_size_total + TS->u_offsets[ns]),
                                     &(sim[ns].solver),
                                     &(sim[ns].mpi),
                                     stagetime);
        }
      }

      for (ns = 0; ns < nsims; ns++) {
          TS->RHSFunction( (TS->gpu_Udot + stage*TS->u_size_total + TS->u_offsets[ns]),
                           (TS->gpu_U + stage*TS->u_size_total + TS->u_offsets[ns]),
                           &(sim[ns].solver),
                           &(sim[ns].mpi),
                           stagetime  );
      }

      gpuArraySetValue(TS->gpu_BoundaryFlux + stage*TS->bf_size_total, TS->bf_size_total, 0.0);

      for (ns = 0; ns < nsims; ns++) {
          gpuArrayCopy1D(  sim[ns].solver.StageBoundaryIntegral,
                           (TS->gpu_BoundaryFlux + stage*TS->bf_size_total + TS->bf_offsets[ns]),
                            TS->bf_sizes[ns]  );

      }

    }

    /* Step completion */
    for (stage = 0; stage < params->nstages; stage++) {

      for (ns = 0; ns < nsims; ns++) {
          gpuArrayAXPY(  (TS->gpu_Udot + stage*TS->u_size_total + TS->u_offsets[ns]),
                          (TS->dt * params->b[stage]),
                          (sim[ns].solver.gpu_u),
                          (TS->u_sizes[ns]) );
          gpuArrayAXPY(  (TS->gpu_BoundaryFlux + stage*TS->bf_size_total + TS->bf_offsets[ns]),
                          (TS->dt * params->b[stage]),
                          (sim[ns].solver.StepBoundaryIntegral),
                          (TS->bf_sizes[ns]) );

      }

    }

  } else {
#endif

    /* Calculate stage values */
    for (stage = 0; stage < params->nstages; stage++) {

      double stagetime = TS->waqt + params->c[stage]*TS->dt;

      for (ns = 0; ns < nsims; ns++) {
        _ArrayCopy1D_(  sim[ns].solver.u,
                        (TS->U[stage] + TS->u_offsets[ns]),
                        (TS->u_sizes[ns]) );
      }

      for (i = 0; i < stage; i++) {
        _ArrayAXPY_(  TS->Udot[i],
                      (TS->dt * params->A[stage*params->nstages+i]),
                      TS->U[stage],
                      TS->u_size_total );
      }

      for (ns = 0; ns < nsims; ns++) {
        if (sim[ns].solver.PreStage) {
          fprintf(stderr,"ERROR in TimeRK(): Call to solver->PreStage() commented out!\n");
          return 1;
  //        sim[ns].solver.PreStage( stage,
  //                                 (TS->U),
  //                                 &(sim[ns].solver),
  //                                 &(sim[ns].mpi),
  //                                 stagetime ); CHECKERR(ierr);
        }
      }

      for (ns = 0; ns < nsims; ns++) {
        if (sim[ns].solver.PostStage) {
          sim[ns].solver.PostStage(  (TS->U[stage] + TS->u_offsets[ns]),
                                     &(sim[ns].solver),
                                     &(sim[ns].mpi),
                                     stagetime); CHECKERR(ierr);
        }
      }

      for (ns = 0; ns < nsims; ns++) {
        TS->RHSFunction( (TS->Udot[stage] + TS->u_offsets[ns]),
                         (TS->U[stage] + TS->u_offsets[ns]),
                         &(sim[ns].solver),
                         &(sim[ns].mpi),
                         stagetime);
      }

      _ArraySetValue_(TS->BoundaryFlux[stage], TS->bf_size_total, 0.0);
      for (ns = 0; ns < nsims; ns++) {
        _ArrayCopy1D_(  sim[ns].solver.StageBoundaryIntegral,
                        (TS->BoundaryFlux[stage] + TS->bf_offsets[ns]),
                        TS->bf_sizes[ns] );
      }

    }

    /* Step completion */
    for (stage = 0; stage < params->nstages; stage++) {

      for (ns = 0; ns < nsims; ns++) {
        _ArrayAXPY_(  (TS->Udot[stage] + TS->u_offsets[ns]),
                      (TS->dt * params->b[stage]),
                      (sim[ns].solver.u),
                      (TS->u_sizes[ns]) );
        _ArrayAXPY_(  (TS->BoundaryFlux[stage] + TS->bf_offsets[ns]),
                      (TS->dt * params->b[stage]),
                      (sim[ns].solver.StepBoundaryIntegral),
                      (TS->bf_sizes[ns]) );
      }

    }

#if defined(HAVE_CUDA)
  }
#endif

  return 0;
}

