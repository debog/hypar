/*! @file main.c
 *  @brief Main driver.
 * The main driver function that calls the initialization, solving, and cleaning up functions.
 *  @author Debojyoti Ghosh
*/

/*! @mainpage

* HyPar - Hyperbolic-Parabolic Partial Differential Equations Solver: 
* A finite-difference algorithm to solve hyperbolic-parabolic equations 
* (with source term). The hyperbolic terms are discretized using a conservative 
* finite-difference scheme (eg: 1st order UPWIND, 3rd order MUSCL, 5th order WENO, 
* 5th order CRWENO). The parabolic terms are discretized either using a conservative 
* or a non-conservative scheme. Time integration is carried out using the PETSc TS 
* library. If compiled without PETSc, the first order Euler and some higher order 
* multi-stage Runge-Kutta schemes are available. Examples of physical models include 
* the linear advection-diffusion-reaction, Euler and Navier-Stokes equations, 
* Fokker-Planck equations for power systems, etc. The code can be compiled in serial 
* as well as in parallel (MPI). For more details, see README.

*  @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#ifdef with_petsc
#include <petscinterface.h>
#endif
#include <mpivars.h>
#include <hypar.h>

static const char help[] = "HyPar - A finite-difference algorithm for solving hyperbolic-parabolic PDEs";

/*!
 * \brief Main driver
 *
 * The main driver function that calls the initialization, solving, and cleaning up functions.
*/
int main(int argc,char **argv)
{
  MPIVariables    mpi;
  HyPar           solver;
  int             ierr = 0, d;
  struct timeval  main_start, solve_start;
  struct timeval  main_end  , solve_end  ;
#ifdef with_petsc
  PetscBool       use_petscts;
#endif

#ifdef serial
  mpi.rank  = 0;
  mpi.nproc = 1;
  mpi.world = 0;
  mpi.comm  = NULL;
  printf("HyPar - Serial Version\n");
#else
  MPI_Init(&argc,&argv);
  MPI_Comm_dup (MPI_COMM_WORLD,&mpi.world);
  MPI_Comm_rank(mpi.world,&mpi.rank );
  MPI_Comm_size(mpi.world,&mpi.nproc);
  if (!mpi.rank) printf("HyPar - Parallel (MPI) version with %d processes\n",mpi.nproc);
#endif

#ifdef with_petsc
  PetscInitialize(&argc,&argv,(char*)0,help);
  if (!mpi.rank) printf("Compiled with PETSc time integration.\n");

  use_petscts = PETSC_FALSE; /* default value */
  ierr = PetscOptionsGetBool(PETSC_NULL,"-use-petscts" ,&use_petscts ,PETSC_NULL); CHKERRQ(ierr);
  solver.use_petscTS  = use_petscts;
#endif

  gettimeofday(&main_start,NULL);

  /* Read Inputs */
  ierr = ReadInputs(&solver,&mpi);
  if (ierr) {
    printf("Error: ReadInputs() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* Initialize and allocate arrays */
  ierr = Initialize(&solver,&mpi);
  if (ierr) {
    printf("Error: Initialize() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* read and set grid & initial solution */
  ierr = InitialSolution(&solver,&mpi);
  if (ierr) {
    printf("Error: InitialSolution() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* Initialize domain boundaries */
  ierr = InitializeBoundaries(&solver,&mpi);
  if (ierr) {
    printf("Error: InitializeBoundaries() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* Initialize solvers */
  ierr = InitializeSolvers(&solver,&mpi);
  if (ierr) {
    printf("Error: InitializeSolvers() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* Initialize physics */
  ierr = InitializePhysics(&solver,&mpi);
  if (ierr) {
    printf("Error: InitializePhysics() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* Initializations complete */
  
  /* Run the solver */
#ifndef serial
  MPI_Barrier(mpi.world);
#endif
  gettimeofday(&solve_start,NULL);
#ifdef with_petsc
  if (solver.use_petscTS == PETSC_TRUE) {
    /* Use PETSc time-integration */
    ierr = SolvePETSc(&solver,&mpi);
    if (ierr) {
      printf("Error: SolvePETSc() returned with status %d on process %d.\n",ierr,mpi.rank);
      return(ierr);
    }
  } else {
    /* Use native time-integration */
    ierr = Solve(&solver,&mpi);
    if (ierr) {
      printf("Error: Solve() returned with status %d on process %d.\n",ierr,mpi.rank);
      return(ierr);
    }
  }
#else 
  /* Use native time-integration */
  ierr = Solve(&solver,&mpi);
  if (ierr) {
    printf("Error: Solve() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
#endif
  gettimeofday(&solve_end,NULL);
#ifndef serial
  MPI_Barrier(mpi.world);
#endif
  gettimeofday(&main_end,NULL);

  /* calculate solver and total runtimes */
  long long walltime;
  walltime = (  (main_end.tv_sec * 1000000   + main_end.tv_usec  ) 
              - (main_start.tv_sec * 1000000 + main_start.tv_usec));
  double main_runtime = (double) walltime / 1000000.0;
  ierr = MPIMax_double(&main_runtime,&main_runtime,1,&mpi.world); if(ierr) return(ierr);
  walltime = (  (solve_end.tv_sec * 1000000   + solve_end.tv_usec  ) 
              - (solve_start.tv_sec * 1000000 + solve_start.tv_usec));
  double solver_runtime = (double) walltime / 1000000.0;
  ierr = MPIMax_double(&solver_runtime,&solver_runtime,1,&mpi.world); if(ierr) return(ierr);

  if (!mpi.rank) {
    FILE *out; 
    /* write out solution errors and wall times to file */
    out = fopen("errors.dat","w");
    for (d=0; d<solver.ndims; d++) fprintf(out,"%4d ",solver.dim_global[d]);
    for (d=0; d<solver.ndims; d++) fprintf(out,"%4d ",mpi.iproc[d]);
    fprintf(out,"%1.16E  ",solver.dt);
    fprintf(out,"%1.16E %1.16E %1.16E   ",solver.error[0],solver.error[1],solver.error[2]);
    fprintf(out,"%1.16E %1.16E\n",solver_runtime,main_runtime);
    fclose(out);
    /* write out conservation errors to file */
    out = fopen("conservation.dat","w");
    for (d=0; d<solver.ndims; d++) fprintf(out,"%4d ",solver.dim_global[d]);
    for (d=0; d<solver.ndims; d++) fprintf(out,"%4d ",mpi.iproc[d]);
    fprintf(out,"%1.16E  ",solver.dt);
    for (d=0; d<solver.nvars; d++) fprintf(out,"%1.16E ",solver.ConservationError[d]);
    fprintf(out,"\n");
    fclose(out);
    /* write out function call counts to file */
    out = fopen("function_counts.dat","w");
    fprintf(out,"%d\n",solver.n_iter);
    fprintf(out,"%d\n",solver.count_hyp);
    fprintf(out,"%d\n",solver.count_par);
    fprintf(out,"%d\n",solver.count_sou);
#ifdef with_petsc
    fprintf(out,"%d\n",solver.count_RHSFunction);
    fprintf(out,"%d\n",solver.count_IFunction);
    fprintf(out,"%d\n",solver.count_IJacobian);
    fprintf(out,"%d\n",solver.count_IJacFunction);
#endif
    fclose(out);
    /* print solution errors, conservation errors, and wall times to screen */
    printf("Computed errors:\n");
    printf("  L1         Error           : %1.16E\n",solver.error[0]);
    printf("  L2         Error           : %1.16E\n",solver.error[1]);
    printf("  Linfinity  Error           : %1.16E\n",solver.error[2]);
    printf("Conservation Errors:\n");
    for (d=0; d<solver.nvars; d++) printf("\t%1.16E\n",solver.ConservationError[d]);
    printf("Solver runtime (in seconds): %1.16E\n",solver_runtime);
    printf("Total  runtime (in seconds): %1.16E\n",main_runtime);
  }

  /* Cleaning up */
  ierr = Cleanup(&solver,&mpi);
  if (ierr) {
    printf("Error: CleanUp() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  if (!mpi.rank) printf("Finished.\n");

#ifdef with_petsc
  PetscFinalize();
#endif

#ifndef serial
  MPI_Comm_free(&mpi.world);
  MPI_Finalize();
#endif
  return(0);
}
