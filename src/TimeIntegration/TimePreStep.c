/*! @file TimePreStep.c
    @brief Pre-time-step function
    @author Debojyoti Ghosh, Youngdae Kim
*/

#include <basic.h>
#if defined(HAVE_CUDA)
#include <arrayfunctions_gpu.h>
#else
#include <arrayfunctions.h>
#endif
#include <timeintegration.h>
#include <mpivars.h>
#include <simulation_object.h>

/*!
  Pre-time-step function: This function is called before each time
  step. Some notable things this does are:
  + Computes CFL and diffusion numbers.
  + Call the physics-specific pre-time-step function, if defined.
*/
int TimePreStep(void *ts /*!< Object of type #TimeIntegration */ )
{
  TimeIntegration*  TS  = (TimeIntegration*) ts;
  _DECLARE_IERR_;

  SimulationObject* sim = (SimulationObject*) TS->simulation;
  int ns, nsims = TS->nsims;

  gettimeofday(&TS->iter_start_time,NULL);

  for (ns = 0; ns < nsims; ns++) {

    HyPar*        solver = &(sim[ns].solver);
    MPIVariables* mpi    = &(sim[ns].mpi);

    double *u = NULL;
#if defined(HAVE_CUDA)
    if (solver->use_gpu) {
      u = solver->gpu_u;
    } else{
#endif
      u = solver->u;
#if defined(HAVE_CUDA)
    }
#endif

    /* apply boundary conditions and exchange data over MPI interfaces */

    solver->ApplyBoundaryConditions( solver,
                                     mpi,
                                     u,
                                     NULL,
                                     TS->waqt );

    solver->ApplyIBConditions( solver,
                               mpi,
                               u,
                               TS->waqt );

#if defined(HAVE_CUDA)
    if (solver->use_gpu) {
      gpuMPIExchangeBoundariesnD( solver->ndims,
                                  solver->nvars,
                                  solver->gpu_dim_local,
                                  solver->ghosts,
                                  mpi,
                                  u );
    } else {
#endif
      MPIExchangeBoundariesnD( solver->ndims,
                               solver->nvars,
                               solver->dim_local,
                               solver->ghosts,
                               mpi,
                               u );
#if defined(HAVE_CUDA)
    }
#endif

    if ((TS->iter+1)%solver->screen_op_iter == 0) {

#if defined(HAVE_CUDA)
      if (!solver->use_gpu) {
#endif
        _ArrayCopy1D_(  solver->u,
                        (TS->u + TS->u_offsets[ns]),
                        (solver->npoints_local_wghosts*solver->nvars) );
#if defined(HAVE_CUDA)
      }
#endif

      /* compute max CFL and diffusion number over the domain */
      if (solver->ComputeCFL) {
        double local_max_cfl  = -1.0;
        local_max_cfl  = solver->ComputeCFL (solver,mpi,TS->dt,TS->waqt);
        MPIMax_double(&TS->max_cfl ,&local_max_cfl ,1,&mpi->world);
      } else {
        TS->max_cfl = -1;
      }
      if (solver->ComputeDiffNumber) {
        double local_max_diff = -1.0;
        local_max_diff = solver->ComputeDiffNumber (solver,mpi,TS->dt,TS->waqt);
        MPIMax_double(&TS->max_diff,&local_max_diff,1,&mpi->world);
      } else {
        TS->max_diff = -1;
      }

    }

    /* set the step boundary flux integral value to zero */
#if defined(HAVE_CUDA)
    if (solver->use_gpu) {
      gpuArraySetValue(solver->StepBoundaryIntegral,2*solver->ndims*solver->nvars,0.0);
    } else {
#endif
      _ArraySetValue_(solver->StepBoundaryIntegral,2*solver->ndims*solver->nvars,0.0);
#if defined(HAVE_CUDA)
    }
#endif

    if (solver->PreStep) {
      solver->PreStep(u,solver,mpi,TS->waqt);
    }

  }

  return 0;
}
