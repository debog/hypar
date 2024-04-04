/*! @file PetscPostTimeStep.cpp
    @brief Post-time-step function
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdio.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars_cpp.h>
#include <simulation_object.h>
#include <petscinterface.h>

int OutputSolution (void*,int,double);   /*!< Write solutions to file */

#undef __FUNCT__
#define __FUNCT__ "PetscPostTimeStep"

/*! Function called after every time step */
PetscErrorCode PetscPostTimeStep(TS ts /*!< Time integrator object */)
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

  Vec Y;
  TSGetSolution(ts,&Y);

  double waqt;
  TSGetTime(ts,&waqt);
  int iter;
  TSGetStepNumber(ts,&iter);
  double dt;
  TSGetTimeStep(ts,&dt);

  context->tic++;

  double max_cfl = 0.0;
  double total_norm = 0.0;
  for (int ns = 0; ns < nsims; ns++) {

    HyPar* solver = &(sim[ns].solver);
    MPIVariables* mpi = &(sim[ns].mpi);

    /* get solution */
    TransferVecFromPETSc(solver->u,Y,context,ns,context->offsets[ns]);

    /* Call any physics-specific post-step function */
    if (solver->PostStep) solver->PostStep(solver->u,solver,mpi,waqt,iter);

    /* Calculate CFL and diffusion number */
    double local_max_cfl  = -1.0, global_max_cfl  = -1.0;
    if (solver->ComputeCFL) local_max_cfl = solver->ComputeCFL(solver,mpi,dt,waqt);
    MPIMax_double(&global_max_cfl ,&local_max_cfl ,1,&mpi->world);
    if (global_max_cfl > max_cfl) max_cfl = global_max_cfl;

    /* Calculate norm of the change in the solution for this time step */
    _ArrayAXPY_(solver->u,-1.0,solver->u0,(solver->npoints_local_wghosts*solver->nvars));
    double sum = ArraySumSquarenD(solver->nvars,solver->ndims,solver->dim_local,
                                  solver->ghosts,solver->index,solver->u0);
    double global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
    double norm = sqrt((global_sum/(double)solver->npoints_global));
    total_norm += norm;

    if (!strcmp(solver->ConservationCheck,"yes")) {
      /* calculate volume integral of the solution at this time step */
      solver->VolumeIntegralFunction(solver->VolumeIntegral,solver->u,solver,mpi);
      /* calculate surface integral of the flux at this time step */
      solver->BoundaryIntegralFunction(solver,mpi);
      /* calculate the conservation error at this time step       */
      solver->CalculateConservationError(solver,mpi);
    }

  }

  gettimeofday(&context->iter_end_time,NULL);
  long long walltime;
  walltime = (  (context->iter_end_time.tv_sec * 1000000 + context->iter_end_time.tv_usec)
              - (context->iter_start_time.tv_sec * 1000000 + context->iter_start_time.tv_usec));
  context->iter_wctime = (double) walltime / 1000000.0;
  context->ti_runtime += context->iter_wctime;

  if ((!context->rank) && (iter%sim[0].solver.screen_op_iter == 0)) {

    if (nsims > 1) {
      printf("--\n");
      printf("iter=%7d,  t=%1.3e\n", iter, waqt);
      printf("  dt=%1.3E  ", dt);
      printf("  CFL=%1.3E, ", max_cfl);
      printf("  norm=%1.4E\n", total_norm);
      printf("  wctime=%1.1E (s)\n",context->iter_wctime);
    } else {
      printf("iter=%7d  ", iter);
      printf("dt=%1.3E  ", dt);
      printf("t=%1.3E  ", waqt);
      printf("CFL=%1.3E  ", max_cfl );
      printf("norm=%1.4E  ", total_norm);
      printf("wctime: %1.1E (s)  ",context->iter_wctime);
    }

    /* calculate and print conservation error */
    if (!strcmp(sim[0].solver.ConservationCheck,"yes")) {

      double error = 0;
      for  (int ns = 0; ns < nsims; ns++) {
        for (int v=0; v<sim[ns].solver.nvars; v++) {
          error += (sim[ns].solver.ConservationError[v] * sim[ns].solver.ConservationError[v]);
        }
      }
      error = sqrt(error);
      printf("  cons_err=%1.4E\n", error);

    } else {

      if (nsims == 1) printf("\n");

    }

    /* print physics-specific info, if available */
    for( int ns = 0; ns < nsims; ns++) {
      if (sim[ns].solver.PrintStep) {
        if (nsims > 1) printf("Physics-specific output for domain %d:\n", ns);
        sim[ns].solver.PrintStep( &(sim[ns].solver),
                                  &(sim[ns].mpi),
                                  waqt );
      }
    }

    if (nsims > 1) {
      printf("--\n");
      printf("\n");
    }

  }

  /* Write intermediate solution to file */
  if (iter%sim[0].solver.file_op_iter == 0) {
    for (int ns = 0; ns < nsims; ns++) {
      HyPar* solver = &(sim[ns].solver);
      MPIVariables* mpi = &(sim[ns].mpi);
      if (solver->PhysicsOutput) solver->PhysicsOutput(solver,mpi,waqt);
    }
    OutputSolution(sim,nsims,waqt);
#ifdef with_librom
    context->op_times_arr.push_back(waqt);
#endif
    context->tic=0;
  }

  PetscFunctionReturn(0);
}

#endif
