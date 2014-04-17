#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#ifdef with_petsc
#include <petscinterface.h>
#endif
#include <mpivars.h>
#include <hypar.h>

static const char help[] = "HyPar - A finite-difference algorithm for solving hyperbolic-parabolic PDEs";

int main(int argc,char **argv)
{
  MPIVariables    mpi;
  HyPar           solver;
  int             ierr = 0,d;
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
  if (!mpi.rank) printf("compiled with PETSc time integration.\n");

  use_petscts = PETSC_FALSE; /* default value */
  ierr = PetscOptionsGetBool(PETSC_NULL,"-use-petscts" ,&use_petscts ,PETSC_NULL); CHKERRQ(ierr);
  solver.use_petscTS  = use_petscts;
#else
  printf(".\n");
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
  
  /* Write an initial solution file */
  ierr = OutputSolution(&solver,&mpi);
  if (ierr) {
    printf("Error: OutputSolution() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }

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

  /* Write final solution file */
  ierr = OutputSolution(&solver,&mpi);
  if (ierr) {
    printf("Error: OutputSolution() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }

  /* Calculate error if exact solution is available */
  ierr = CalculateError(&solver,&mpi);
  if (ierr) {
    printf("Error: CalculateError() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }

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

  /* print error and walltime to file and on screen */
  if (!mpi.rank) {
    FILE *out; out = fopen("errors.dat","w");
    for (d=0; d<solver.ndims; d++) fprintf(out,"%4d ",solver.dim_global[d]);
    for (d=0; d<solver.ndims; d++) fprintf(out,"%4d ",mpi.iproc[d]);
    fprintf(out,"%1.16E %1.16E %1.16E   ",solver.error[0],solver.error[1],solver.error[2]);
    fprintf(out,"%1.16E %1.16E\n",solver_runtime,main_runtime);
    fclose(out);
    printf("L1         Error           : %1.16E\n",solver.error[0]);
    printf("L2         Error           : %1.16E\n",solver.error[1]);
    printf("Linfinity  Error           : %1.16E\n",solver.error[2]);
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


#ifndef serial
  MPI_Comm_free(&mpi.world);
  MPI_Finalize();
#endif
  return(0);
}
