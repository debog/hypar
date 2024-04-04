/*! @file PetscPostStage.cpp
    @author Debojyoti Ghosh
    @brief Post-time-integration-stage function.
*/
#ifdef with_petsc

#include <stdio.h>
#include <basic.h>
#include <simulation_object.h>
#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscPostStage"

/*! Function called after every stage in a multi-stage time-integration method */
PetscErrorCode PetscPostStage(  TS        ts,         /*!< Time integrator of PETSc type TS */
                                PetscReal stagetime,  /*!< Current stage time */
                                PetscInt  stageindex, /*!< Stage */
                                Vec       *Y          /*!< Stage solutions (all stages) -
                                                           be careful what you access */ )
{
  PETScContext* context(nullptr);

  PetscFunctionBegin;

  TSGetApplicationContext(ts,&context);
  if (!context) {
    fprintf(stderr,"Error in PetscPreTimeStep: Null context!\n");
    return(1);
  }
  SimulationObject* sim = (SimulationObject*) context->simobj;
  int nsims = context->nsims;

  TSType time_scheme;
  TSGetType(ts,&time_scheme);

  TSGetTimeStep(ts,&(context->dt));
  context->stage_index = stageindex;
  if (context->stage_times.size() == stageindex) {
    context->stage_times.push_back(stagetime/context->dt);
  }

  for (int ns = 0; ns < nsims; ns++) {

    HyPar* solver = &(sim[ns].solver);
    MPIVariables* mpi = &(sim[ns].mpi);

    /* get solution */
    TransferVecFromPETSc(solver->u,Y[stageindex],context,ns,context->offsets[ns]);

    /* apply immersed boundaries */
    solver->ApplyIBConditions(solver,mpi,solver->u,stagetime);

    /* If using a non-linear scheme with ARKIMEX methods,
       compute the non-linear finite-difference operator */
    if (!strcmp(time_scheme,TSARKIMEX)) {
      solver->NonlinearInterp(  solver->u,
                                solver,
                                mpi,
                                (double)stagetime,
                                solver->FFunction );
    }

    /* Call any physics-specific post-stage function, if available */
    if (solver->PostStage) {
      solver->PostStage(solver->u,solver,mpi,stagetime);
    }

    TransferVecToPETSc(solver->u,Y[stageindex],context,ns,context->offsets[ns]);

  }

  PetscFunctionReturn(0);
}

#endif
