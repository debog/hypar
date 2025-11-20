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
PetscErrorCode PetscPostStage(  TS        a_ts,         /*!< Time integrator of PETSc type TS */
                                PetscReal a_stagetime,  /*!< Current stage time */
                                PetscInt  a_stageindex, /*!< Stage */
                                Vec       *a_Y          /*!< Stage solutions (all stages) -
                                                           be careful what you access */ )
{
  PETScContext* context(nullptr);

  PetscFunctionBegin;

  TSGetApplicationContext(a_ts,&context);
  if (!context) {
    fprintf(stderr,"Error in PetscPreTimeStep: Null context!\n");
    return(1);
  }
  SimulationObject* sim = (SimulationObject*) context->m_simobj;
  int nsims = context->m_nsims;

  TSType time_scheme;
  TSGetType(a_ts,&time_scheme);

  TSGetTimeStep(a_ts,&(context->m_dt));
  context->m_stage_index = a_stageindex;
  if (context->m_stage_times.size() == a_stageindex) {
    context->m_stage_times.push_back(a_stagetime/context->m_dt);
  }

  for (int ns = 0; ns < nsims; ns++) {

    HyPar* solver = &(sim[ns].solver);
    MPIVariables* mpi = &(sim[ns].mpi);

    /* get solution */
    TransferVecFromPETSc(solver->m_u,a_Y[a_stageindex],context,ns,context->m_offsets[ns]);

    /* apply immersed boundaries */
    solver->ApplyIBConditions(solver,mpi,solver->m_u,a_stagetime);

    /* If using a non-linear scheme with ARKIMEX methods,
       compute the non-linear finite-difference operator */
    if (!strcmp(time_scheme,TSARKIMEX)) {
      solver->NonlinearInterp(  solver->m_u,
                                solver,
                                mpi,
                                (double)a_stagetime,
                                solver->FFunction );
    }

    /* Call any physics-specific post-stage function, if available */
    if (solver->PostStage) {
      solver->PostStage(solver->m_u,solver,mpi,a_stagetime);
    }

    TransferVecToPETSc(solver->m_u,a_Y[a_stageindex],context,ns,context->m_offsets[ns]);

  }

  PetscFunctionReturn(0);
}

#endif
