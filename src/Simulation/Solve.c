/*! @file Solve.c
    @author Debojyoti Ghosh
    @brief  Solve the governing equations in time
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <common.h>
#include <io.h>
#include <timeintegration.h>
#include <simulation.h>

#ifdef compute_rhs_operators
int ComputeRHSOperators(void*,void*,double);
#endif

/*! This function integrates the semi-discrete ODE (obtained from discretizing the 
    PDE in space) using natively implemented time integration methods. It initializes 
    the time integration object, iterates the simulation for the required number of 
    time steps, and calculates the errors. After the specified number of iterations, 
    it writes out some information to the screen and the solution to a file.
*/
int Solve(  void  *s,     /*!< Array of simulation objects of type #SimulationObject */
            int   nsims,  /*!< number of simulation objects */
            int   rank,   /*!< MPI rank of this process */
            int   nproc   /*!< Number of MPI processes */
         )
{
  SimulationObject* sim = (SimulationObject*) s;
  int ns;
  int tic     = 0;
  _DECLARE_IERR_;

  /* write out iblank to file for visualization */
  for (ns = 0; ns < nsims; ns++) {
    if (sim[ns].solver.flag_ib) {

      char fname_root[_MAX_STRING_SIZE_] = "iblank";
      if (nsims > 1) {
        char index[_MAX_STRING_SIZE_];
        GetStringFromInteger(ns, index, (int)log10((nsims)+1));
        strcat(fname_root, "_");
        strcat(fname_root, index);
      }

      WriteArray( sim[ns].solver.ndims,
                  1,
                  sim[ns].solver.dim_global,
                  sim[ns].solver.dim_local,
                  sim[ns].solver.ghosts,
                  sim[ns].solver.x,
                  sim[ns].solver.iblank,
                  &(sim[ns].solver),
                  &(sim[ns].mpi),
                  fname_root );
    }
  }

  /* Define and initialize the time-integration object */
  TimeIntegration TS;
  if (!rank) printf("Setting up time integration.\n");
  IERR TimeInitialize(sim, nsims, rank, nproc, &TS); CHECKERR(ierr);

  if (!rank) printf("Solving in time (from %d to %d iterations)\n",TS.restart_iter,TS.n_iter);
  for (TS.iter = TS.restart_iter; TS.iter < TS.n_iter; TS.iter++) {

    /* Write initial solution to file if this is the first iteration */
    if (!TS.iter) { 
      for (ns = 0; ns < nsims; ns++) {
        if (sim[ns].solver.PhysicsOutput) {
          sim[ns].solver.PhysicsOutput( &(sim[ns].solver),
                                        &(sim[ns].mpi) );
        }
      }
      IERR OutputSolution(sim, nsims); CHECKERR(ierr); 
    }

    /* Call pre-step function */
    IERR TimePreStep  (&TS); CHECKERR(ierr);
#ifdef compute_rhs_operators
    /* compute and write (to file) matrix operators representing the right-hand side */
//    if (((TS.iter+1)%solver->file_op_iter == 0) || (!TS.iter)) 
//      { IERR ComputeRHSOperators(solver,mpi,TS.waqt); CHECKERR(ierr); }
#endif

    /* Step in time */
    IERR TimeStep     (&TS); CHECKERR(ierr);

    /* Call post-step function */
    IERR TimePostStep (&TS); CHECKERR(ierr);

    /* Print information to screen */
    IERR TimePrintStep(&TS); CHECKERR(ierr);
    tic++;

    /* Write intermediate solution to file */
    if ((TS.iter+1)%sim[0].solver.file_op_iter == 0) { 
      for (ns = 0; ns < nsims; ns++) {
        if (sim[ns].solver.PhysicsOutput) {
          sim[ns].solver.PhysicsOutput( &(sim[ns].solver),
                                        &(sim[ns].mpi) );
        }
      }
      IERR OutputSolution(sim, nsims); CHECKERR(ierr); 
      tic = 0; 
    }

  }

  /* write a final solution file, if last iteration did not write one */
  if (tic || (!TS.n_iter)) { 
    for (ns = 0; ns < nsims; ns++) {
      if (sim[ns].solver.PhysicsOutput) {
        sim[ns].solver.PhysicsOutput( &(sim[ns].solver),
                                      &(sim[ns].mpi) );
      }
    }
    IERR OutputSolution(sim, nsims); CHECKERR(ierr); 
  }

  if (!rank) {
    printf("Completed time integration (Final time: %f).\n",TS.waqt);
    if (nsims > 1) printf("\n");
  }

  /* calculate error if exact solution has been provided */
  for (ns = 0; ns < nsims; ns++) {
    IERR CalculateError(&(sim[ns].solver),
                        &(sim[ns].mpi) ); CHECKERR(ierr);
  }
  IERR TimeCleanup(&TS); CHECKERR(ierr);

  return(0);
}
