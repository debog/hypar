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
int TimePreStep(void *a_ts /*!< Object of type #TimeIntegration */ )
{
  TimeIntegration*  TS  = (TimeIntegration*) a_ts;
  _DECLARE_IERR_;

  SimulationObject* sim = (SimulationObject*) TS->m_simulation;
  int ns, nsims = TS->m_nsims;

  gettimeofday(&TS->m_iter_start_time,NULL);

  for (ns = 0; ns < nsims; ns++) {

    HyPar*        solver = &(sim[ns].solver);
    MPIVariables* mpi    = &(sim[ns].mpi);

    double *u = NULL;
#if defined(HAVE_CUDA)
    if (solver->m_use_gpu) {
      u = solver->m_gpu_u;
    } else{
#endif
      u = solver->m_u;
#if defined(HAVE_CUDA)
    }
#endif

    /* apply boundary conditions and exchange data over MPI interfaces */

    solver->ApplyBoundaryConditions( solver,
                                     mpi,
                                     u,
                                     NULL,
                                     TS->m_waqt );

    solver->ApplyIBConditions( solver,
                               mpi,
                               u,
                               TS->m_waqt );

#if defined(HAVE_CUDA)
    if (solver->m_use_gpu) {
      gpuMPIExchangeBoundariesnD( solver->m_ndims,
                                  solver->m_nvars,
                                  solver->m_gpu_dim_local,
                                  solver->m_ghosts,
                                  mpi,
                                  u );
    } else {
#endif
      MPIExchangeBoundariesnD( solver->m_ndims,
                               solver->m_nvars,
                               solver->m_dim_local,
                               solver->m_ghosts,
                               mpi,
                               u );
#if defined(HAVE_CUDA)
    }
#endif

    if ((TS->m_iter+1)%solver->m_screen_op_iter == 0) {

#if defined(HAVE_CUDA)
      if (!solver->m_use_gpu) {
#endif
        _ArrayCopy1D_(  solver->m_u,
                        (TS->m_u + TS->m_u_offsets[ns]),
                        (solver->m_npoints_local_wghosts*solver->m_nvars) );
#if defined(HAVE_CUDA)
      }
#endif

      /* compute max CFL and diffusion number over the domain */
      if (solver->ComputeCFL) {
        double local_max_cfl  = -1.0;
        local_max_cfl  = solver->ComputeCFL (solver,mpi,TS->m_dt,TS->m_waqt);
        MPIMax_double(&TS->m_max_cfl ,&local_max_cfl ,1,&mpi->m_world);
      } else {
        TS->m_max_cfl = -1;
      }
      if (solver->ComputeDiffNumber) {
        double local_max_diff = -1.0;
        local_max_diff = solver->ComputeDiffNumber (solver,mpi,TS->m_dt,TS->m_waqt);
        MPIMax_double(&TS->m_max_diff,&local_max_diff,1,&mpi->m_world);
      } else {
        TS->m_max_diff = -1;
      }

    }

    /* set the step boundary flux integral value to zero */
#if defined(HAVE_CUDA)
    if (solver->m_use_gpu) {
      gpuArraySetValue(solver->m_step_boundary_integral,2*solver->m_ndims*solver->m_nvars,0.0);
    } else {
#endif
      _ArraySetValue_(solver->m_step_boundary_integral,2*solver->m_ndims*solver->m_nvars,0.0);
#if defined(HAVE_CUDA)
    }
#endif

    if (solver->PreStep) {
      solver->PreStep(u,solver,mpi,TS->m_waqt);
    }

  }

  return 0;
}
