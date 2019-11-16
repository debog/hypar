/*! @file TimePreStep.c
    @brief Pre-time-step function
    @author Debojyoti Ghosh
*/

#include <basic.h>
#include <arrayfunctions.h>
#include <timeintegration.h>
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

  for (ns = 0; ns < nsims; ns++) {

    HyPar*        solver = &(sim[ns].solver);
    MPIVariables* mpi    = &(sim[ns].mpi);
  
    /* apply boundary conditions and exchange data over MPI interfaces */
    IERR solver->ApplyBoundaryConditions( solver,
                                          mpi,
                                          solver->u,
                                          NULL,
                                          TS->waqt ); CHECKERR(ierr);

    IERR solver->ApplyIBConditions( solver,
                                    mpi,
                                    solver->u,
                                    TS->waqt ); CHECKERR(ierr);

    IERR MPIExchangeBoundariesnD( solver->ndims,
                                  solver->nvars,
                                  solver->dim_local,
                                  solver->ghosts,
                                  mpi,
                                  solver->u); CHECKERR(ierr);
  
    if ((TS->iter+1)%solver->screen_op_iter == 0) {

      _ArrayCopy1D_(  solver->u,
                      (TS->u + TS->u_offsets[ns]), 
                      (solver->npoints_local_wghosts*solver->nvars) );
  
      /* compute max CFL and diffusion number over the domain */

      double local_max_cfl  = -1.0;
      double local_max_diff = -1.0;

      if (solver->ComputeCFL       ) local_max_cfl  = solver->ComputeCFL        (solver,mpi,TS->dt,TS->waqt);
      if (solver->ComputeDiffNumber) local_max_diff = solver->ComputeDiffNumber (solver,mpi,TS->dt,TS->waqt);

      IERR MPIMax_double(&TS->max_cfl ,&local_max_cfl ,1,&mpi->world); CHECKERR(ierr);
      IERR MPIMax_double(&TS->max_diff,&local_max_diff,1,&mpi->world); CHECKERR(ierr);

    }
  
    /* set the step boundary flux integral value to zero */
    _ArraySetValue_(solver->StepBoundaryIntegral,2*solver->ndims*solver->nvars,0.0);
  
    if (solver->PreStep)  { IERR solver->PreStep(solver->u,solver,mpi,TS->waqt); CHECKERR(ierr); }

  }

  return(0);
}

