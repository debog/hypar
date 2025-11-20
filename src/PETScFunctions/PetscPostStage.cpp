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
  SimulationObject* sim = (SimulationObject*) context->m_simobj;
  int nsims = context->m_nsims;

  TSType time_scheme;
  TSGetType(ts,&time_scheme);

  TSGetTimeStep(ts,&(context->m_dt));
  context->m_stage_index = stageindex;
  if (context->m_stage_times.size() == stageindex) {
    context->m_stage_times.push_back(stagetime/context->m_dt);
  }

  for (int ns = 0; ns < nsims; ns++) {

    HyPar* solver = &(sim[ns].solver);
    MPIVariables* mpi = &(sim[ns].mpi);

    /* get solution */
    TransferVecFromPETSc(solver->m_u,Y[stageindex],context,ns,context->m_offsets[ns]);

    /* apply immersed boundaries */
    solver->ApplyIBConditions(solver,mpi,solver->m_u,stagetime);

    /* If using a non-linear scheme with ARKIMEX methods,
       compute the non-linear finite-difference operator */
    if (!strcmp(time_scheme,TSARKIMEX)) {
      solver->NonlinearInterp(  solver->m_u,
                                solver,
                                mpi,
                                (double)stagetime,
                                solver->FFunction );
    }

    /* Call any physics-specific post-stage function, if available */
    if (solver->PostStage) {
      solver->PostStage(solver->m_u,solver,mpi,stagetime);
    }

    TransferVecToPETSc(solver->m_u,Y[stageindex],context,ns,context->m_offsets[ns]);

  }

  PetscFunctionReturn(0);
}

#endif
