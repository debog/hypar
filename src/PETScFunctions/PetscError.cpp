#ifdef with_petsc

/*! @file PetscError.cpp
    @author Debojyoti Ghosh
    @brief Compute time integration error estimates, if available
*/

#include <math.h>
#include <string>
#include <arrayfunctions.h>
#include <mpivars_cpp.h>
#include <simulation_object.h>
#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscTimeError"

/*! Compute the norms of the error estimate, if the PETSc time integrator
    has it (for example TSGLEE) */
int PetscTimeError(TS  ts /*!< Time integrator object of PETSc type TS */)
{
  PetscErrorCode ierr;

  PETScContext* context(nullptr);
  ierr = TSGetApplicationContext(ts,&context); CHKERRQ(ierr);
  if (!context) {
    fprintf(stderr,"Error in PetscError: Null context!\n");
    return(1);
  }

  SimulationObject* sim( (SimulationObject*)context->simobj );
  int nsims( context->nsims );

  double dt;
  ierr = TSGetTimeStep(ts,&dt); CHKERRQ(ierr);
  TSType time_scheme;
  ierr = TSGetType(ts,&time_scheme); CHKERRQ(ierr);
  Vec Y;
  ierr = TSGetSolution(ts,&Y); CHKERRQ(ierr);

  for (int ns = 0; ns < nsims; ns++) {
    TransferVecFromPETSc( sim[ns].solver.u,
                          Y,
                          context,
                          ns,
                          context->offsets[ns] );
  }

  if (std::string(time_scheme) == std::string(TSGLEE)) {

    Vec Z;
    ierr = VecDuplicate(Y,&Z); CHKERRQ(ierr);
    ierr = TSGetTimeError(ts,0,&Z);CHKERRQ(ierr);
    for (int ns = 0; ns < nsims; ns++) {
      TransferVecFromPETSc( sim[ns].solver.uref,
                            Z,
                            context,
                            ns,
                            context->offsets[ns] );
    }
    ierr = VecDestroy(&Z); CHKERRQ(ierr);

    for (int ns = 0; ns < nsims; ns++) {

      HyPar* solver = &(sim[ns].solver);
      MPIVariables* mpi = &(sim[ns].mpi);

      int size = solver->npoints_local_wghosts * solver->nvars;
      double  sum = 0.0,
              global_sum = 0.0,
              *Uerr = solver->uref,
              error[3] = {0,0,0};

      /* calculate solution norm for relative errors */
      double sol_norm[3] = {0.0,0.0,0.0};
      /* L1 */
      sum = ArraySumAbsnD   (solver->nvars,solver->ndims,solver->dim_local,
                             solver->ghosts,solver->index,solver->u);
      global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
      sol_norm[0] = global_sum/((double)solver->npoints_global);
      /* L2 */
      sum = ArraySumSquarenD(solver->nvars,solver->ndims,solver->dim_local,
                             solver->ghosts,solver->index,solver->u);
      global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
      sol_norm[1] = sqrt(global_sum/((double)solver->npoints_global));
      /* Linf */
      sum = ArrayMaxnD      (solver->nvars,solver->ndims,solver->dim_local,
                             solver->ghosts,solver->index,solver->u);
      global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->world);
      sol_norm[2] = global_sum;

      /* calculate L1 norm of error */
      sum = ArraySumAbsnD   (solver->nvars,solver->ndims,solver->dim_local,
                             solver->ghosts,solver->index,Uerr);
      global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
      error[0] = global_sum/((double)solver->npoints_global);
      /* calculate L2 norm of error */
      sum = ArraySumSquarenD(solver->nvars,solver->ndims,solver->dim_local,
                             solver->ghosts,solver->index,Uerr);
      global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
      error[1] = sqrt(global_sum/((double)solver->npoints_global));
      /* calculate Linf norm of error */
      sum = ArrayMaxnD      (solver->nvars,solver->ndims,solver->dim_local,
                             solver->ghosts,solver->index,Uerr);
      global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->world);
      error[2] = global_sum;

      if (   (sol_norm[0] > _MACHINE_ZERO_)
          && (sol_norm[1] > _MACHINE_ZERO_)
          && (sol_norm[2] > _MACHINE_ZERO_) ) {
        error[0] /= sol_norm[0];
        error[1] /= sol_norm[1];
        error[2] /= sol_norm[2];
      }

      /* write to file */
      if (!mpi->rank) {
        std::string fname = "glm_err";
        if (nsims > 1) {
          char idx_string[10];
          sprintf(idx_string, "%3d", ns);
          fname += ("_" + std::string(idx_string));
        }
        fname += ".dat";

        FILE *out;
        out = fopen(fname.c_str(),"w");
        fprintf(out,"%1.16E  %1.16E  %1.16E  %1.16E  ",dt,error[0],error[1],error[2]);
        fclose(out);
        if (nsims > 1) {
          printf( "Estimated time integration errors (GLM-GEE time-integration) for simulation domain %d:-\n",
                  ns);
        } else {
          printf("Estimated time integration errors (GLM-GEE time-integration):-\n");
        }
        printf("  L1         Error           : %1.16E\n",error[0]);
        printf("  L2         Error           : %1.16E\n",error[1]);
        printf("  Linfinity  Error           : %1.16E\n",error[2]);
      }
    }
  }

  return 0;
}

#endif
