/*! @file PetscPreTimeStep.cpp
    @brief Pre-time-step function
    @author Debojyoti Ghosh */

#ifdef with_petsc

#include <stdio.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars_cpp.h>
#include <simulation_object.h>
#include <petscinterface.h>

#ifdef with_librom
#include <librom_interface.h>
#endif

int OutputSolution (void*,int,double);   /*!< Write solutions to file */

#undef __FUNCT__
#define __FUNCT__ "PetscPreTimeStep"

/*! Function called before a time step */
PetscErrorCode PetscPreTimeStep(TS ts /*!< Time integration object */)
{
  PETScContext* context(nullptr);

  PetscFunctionBegin;

  TSGetApplicationContext(ts,&context);
  if (!context) {
    fprintf(stderr,"Error in PetscPreTimeStep: Null context!\n");
    return(1);
  }
  gettimeofday(&(context->iter_start_time),NULL);
  SimulationObject* sim = (SimulationObject*) context->simobj;
  int nsims = context->nsims;

  Vec Y;
  TSGetSolution(ts,&Y);

  double waqt;
  TSGetTime(ts,&waqt);
  double dt;
  TSGetTimeStep(ts,&dt);
  int iter;
  TSGetStepNumber(ts,&iter);

  context->dt = dt;
  context->waqt = waqt;
  context->t_start = waqt;

  TSType time_scheme;
  TSGetType(ts,&time_scheme);

  for (int ns = 0; ns < nsims; ns++) {

    HyPar* solver = &(sim[ns].solver);
    MPIVariables* mpi = &(sim[ns].mpi);

    /* get solution */
    TransferVecFromPETSc(solver->u,Y,context,ns,context->offsets[ns]);

    /* save a copy of the solution to compute norm at end of time step */
    _ArrayCopy1D_(solver->u,solver->u0,(solver->npoints_local_wghosts*solver->nvars));

    /* apply boundary conditions and exchange data over MPI interfaces */
    solver->ApplyBoundaryConditions(solver,mpi,solver->u,NULL,waqt);
    solver->ApplyIBConditions(solver,mpi,solver->u,waqt);
    MPIExchangeBoundariesnD(  solver->ndims,
                              solver->nvars,
                              solver->dim_local,
                              solver->ghosts,
                              mpi,
                              solver->u );

    /* Call any physics-specific pre-step function */
    if (solver->PreStep) solver->PreStep(solver->u,solver,mpi,waqt);

    /* If using a non-linear scheme with ARKIMEX methods,
       compute the non-linear finite-difference operator */
    if (!strcmp(time_scheme,TSARKIMEX)) {
      solver->NonlinearInterp(solver->u,solver,mpi,waqt,solver->FFunction);
    }

    /* set the step boundary flux integral value to zero */
    _ArraySetValue_(solver->StepBoundaryIntegral,2*solver->ndims*solver->nvars,0.0);

    TransferVecToPETSc(solver->u,Y,context,ns,context->offsets[ns]);

  }


  if (!iter) {
    for (int ns = 0; ns < nsims; ns++) {
      HyPar* solver = &(sim[ns].solver);
      MPIVariables* mpi = &(sim[ns].mpi);
      if (solver->PhysicsOutput) solver->PhysicsOutput(solver,mpi,waqt);
    }
    OutputSolution(sim, nsims,waqt);
#ifdef with_librom
    context->op_times_arr.push_back(waqt);
#endif
  }

#ifdef with_librom
  if (      (context->rom_mode == _ROM_MODE_TRAIN_)
        &&  (iter%((libROMInterface*)context->rom_interface)->samplingFrequency() == 0)  ) {
    ((libROMInterface*)context->rom_interface)->takeSample( sim, waqt );
  }
#endif

  PetscFunctionReturn(0);
}

#endif
