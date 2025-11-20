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
  by one time step of size #HyPar::m_dt using the forward Euler method
  given by
  \f{align}{
    {\bf U}^{\left(i\right)} &= {\bf u}_n + \Delta t \sum_{j=1}^{i-1} a_{ij} {\bf F}\left({\bf U}^{\left(j\right)}\right), \\
    {\bf u}_{n+1} &= {\bf u}_n + \Delta t \sum_{i=1}^s b_{i} {\bf F}\left({\bf U}^{\left(i\right)}\right),
  \f}
  where the subscript represents the time level, the superscripts represent the stages, \f$\Delta t\f$ is the
  time step size #HyPar::m_dt, and \f${\bf F}\left({\bf u}\right)\f$ is computed by #TimeIntegration::RHSFunction.
  The Butcher tableaux coefficients are \f$a_{ij}\f$ (#ExplicitRKParameters::A) and \f$b_i\f$
  (#ExplicitRKParameters::b).

  Note: In the code #TimeIntegration::Udot is equivalent to \f${\bf F}\left({\bf u}\right)\f$.
*/
int TimeRK(void *ts /*!< Object of type #TimeIntegration */)
{
  TimeIntegration* TS = (TimeIntegration*) ts;
  SimulationObject* sim = (SimulationObject*) TS->m_simulation;
  ExplicitRKParameters *params = (ExplicitRKParameters*) sim[0].solver.m_msti;
  int ns, stage, i, nsims = TS->m_nsims;

#if defined(HAVE_CUDA)
  if (sim[0].solver.m_use_gpu) {

    /* Calculate stage values */
    for (stage = 0; stage < params->nstages; stage++) {

      double stagetime = TS->m_waqt + params->c[stage]*TS->m_dt;

      for (ns = 0; ns < nsims; ns++) {
          gpuArrayCopy1D(  sim[ns].solver.m_gpu_u,
                           (TS->m_gpu_U + stage*TS->m_u_size_total + TS->m_u_offsets[ns]),
                           (TS->m_u_sizes[ns])  );

      }

      for (i = 0; i < stage; i++) {
          gpuArrayAXPY(  TS->m_gpu_Udot + i*TS->m_u_size_total,
                         (TS->m_dt * params->A[stage*params->nstages+i]),
                         TS->m_gpu_U + stage*TS->m_u_size_total,
                         TS->m_u_size_total  );

      }

      for (ns = 0; ns < nsims; ns++) {
        if (sim[ns].solver.PreStage) {
          fprintf(stderr,"ERROR in TimeRK(): Call to solver->PreStage() commented out!\n");
          return 1;
//          sim[ns].solver.PreStage( stage,
//                                   (TS->m_gpu_U),
//                                   &(sim[ns].solver),
//                                   &(sim[ns].mpi),
//                                   stagetime ); CHECKERR(ierr);
        }
      }

      for (ns = 0; ns < nsims; ns++) {
        if (sim[ns].solver.PostStage) {
          sim[ns].solver.PostStage(  (TS->m_gpu_U + stage*TS->m_u_size_total + TS->m_u_offsets[ns]),
                                     &(sim[ns].solver),
                                     &(sim[ns].mpi),
                                     stagetime);
        }
      }

      for (ns = 0; ns < nsims; ns++) {
          TS->RHSFunction( (TS->m_gpu_Udot + stage*TS->m_u_size_total + TS->m_u_offsets[ns]),
                           (TS->m_gpu_U + stage*TS->m_u_size_total + TS->m_u_offsets[ns]),
                           &(sim[ns].solver),
                           &(sim[ns].mpi),
                           stagetime  );
      }

      gpuArraySetValue(TS->m_gpu_BoundaryFlux + stage*TS->m_bf_size_total, TS->m_bf_size_total, 0.0);

      for (ns = 0; ns < nsims; ns++) {
          gpuArrayCopy1D(  sim[ns].solver.m_stage_boundary_integral,
                           (TS->m_gpu_BoundaryFlux + stage*TS->m_bf_size_total + TS->m_bf_offsets[ns]),
                            TS->m_bf_sizes[ns]  );

      }

    }

    /* Step completion */
    for (stage = 0; stage < params->nstages; stage++) {

      for (ns = 0; ns < nsims; ns++) {
          gpuArrayAXPY(  (TS->m_gpu_Udot + stage*TS->m_u_size_total + TS->m_u_offsets[ns]),
                          (TS->m_dt * params->b[stage]),
                          (sim[ns].solver.m_gpu_u),
                          (TS->m_u_sizes[ns]) );
          gpuArrayAXPY(  (TS->m_gpu_BoundaryFlux + stage*TS->m_bf_size_total + TS->m_bf_offsets[ns]),
                          (TS->m_dt * params->b[stage]),
                          (sim[ns].solver.m_step_boundary_integral),
                          (TS->m_bf_sizes[ns]) );

      }

    }

  } else {
#endif

    /* Calculate stage values */
    for (stage = 0; stage < params->nstages; stage++) {

      double stagetime = TS->m_waqt + params->c[stage]*TS->m_dt;

      for (ns = 0; ns < nsims; ns++) {
        _ArrayCopy1D_(  sim[ns].solver.m_u,
                        (TS->m_U[stage] + TS->m_u_offsets[ns]),
                        (TS->m_u_sizes[ns]) );
      }

      for (i = 0; i < stage; i++) {
        _ArrayAXPY_(  TS->m_Udot[i],
                      (TS->m_dt * params->A[stage*params->nstages+i]),
                      TS->m_U[stage],
                      TS->m_u_size_total );
      }

      for (ns = 0; ns < nsims; ns++) {
        if (sim[ns].solver.PreStage) {
          fprintf(stderr,"ERROR in TimeRK(): Call to solver->PreStage() commented out!\n");
          return 1;
  //        sim[ns].solver.PreStage( stage,
  //                                 (TS->m_U),
  //                                 &(sim[ns].solver),
  //                                 &(sim[ns].mpi),
  //                                 stagetime ); CHECKERR(ierr);
        }
      }

      for (ns = 0; ns < nsims; ns++) {
        if (sim[ns].solver.PostStage) {
          sim[ns].solver.PostStage(  (TS->m_U[stage] + TS->m_u_offsets[ns]),
                                     &(sim[ns].solver),
                                     &(sim[ns].mpi),
                                     stagetime); CHECKERR(ierr);
        }
      }

      for (ns = 0; ns < nsims; ns++) {
        TS->RHSFunction( (TS->m_Udot[stage] + TS->m_u_offsets[ns]),
                         (TS->m_U[stage] + TS->m_u_offsets[ns]),
                         &(sim[ns].solver),
                         &(sim[ns].mpi),
                         stagetime);
      }

      _ArraySetValue_(TS->m_BoundaryFlux[stage], TS->m_bf_size_total, 0.0);
      for (ns = 0; ns < nsims; ns++) {
        _ArrayCopy1D_(  sim[ns].solver.m_stage_boundary_integral,
                        (TS->m_BoundaryFlux[stage] + TS->m_bf_offsets[ns]),
                        TS->m_bf_sizes[ns] );
      }

    }

    /* Step completion */
    for (stage = 0; stage < params->nstages; stage++) {

      for (ns = 0; ns < nsims; ns++) {
        _ArrayAXPY_(  (TS->m_Udot[stage] + TS->m_u_offsets[ns]),
                      (TS->m_dt * params->b[stage]),
                      (sim[ns].solver.m_u),
                      (TS->m_u_sizes[ns]) );
        _ArrayAXPY_(  (TS->m_BoundaryFlux[stage] + TS->m_bf_offsets[ns]),
                      (TS->m_dt * params->b[stage]),
                      (sim[ns].solver.m_step_boundary_integral),
                      (TS->m_bf_sizes[ns]) );
      }

    }

#if defined(HAVE_CUDA)
  }
#endif

  return 0;
}

