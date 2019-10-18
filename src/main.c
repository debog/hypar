/*! @file main.c
 *  @brief Main driver.
 * The main driver function that calls the initialization, solving, and cleaning up functions.
 *  @author Debojyoti Ghosh
*/

/*! @mainpage

  @author Debojyoti Ghosh [\b Email: (first name) (dot) (last name) (at) gmail (dot) com, \b Website: http://debog.github.io/]

  HyPar: Hyperbolic-Parabolic (with Source) Partial Differential Equations Solver
  -------------------------------------------------------------------------------

  HyPar is a finite-difference algorithm to solve hyperbolic-parabolic partial differential
  equations (with source terms) on Cartesian grids. It is a unified framework that can handle 
  systems of PDEs with arbitrary number of spatial dimensions and solution components. It 
  provides the spatial discretization and time integration functions, functions to read and 
  write solutions from/to files, as well as functions required to solve the system on parallel 
  (MPI) platforms. The physical models define the physics-specific functions such as the exact 
  forms of the hyperbolic flux, parabolic flux, source terms, upwinding functions, etc.

  Features
  --------
  + Solves <B>hyperbolic-parabolic PDEs with source terms</B>.
  + Allows arbitrary number of <B>spatial dimensions</B> and <B>vector components per grid point</B>.
  + Solves the PDEs over <B>Cartesian</B> grids.
  + Written entirely in C and uses the MPICH library. It also uses OpenMP threads 
    but this is a work-in-progress.
  + Can be <B>compiled with PETSc</B> (http://www.mcs.anl.gov/petsc/), if available, where 
    it can use PETSc's time integration module TS (http://www.mcs.anl.gov/petsc/petsc-current/src/ts/).
  + For 3-dimensional simulations, the <B>immersed boundaries</B> can be used to
    solve over non-Cartesian geometries.

  HyPar has been developed to be scalable, and apart from the usual functionalities to
  solve a system of PDEs on distributed architectures, it provides scalable file I/O
  functions. It has been tested on several platforms, including DOE Leadership-class
  supercomputers, with up to ~0.5 million MPI ranks.

  Download
  --------
  The code is available at: https://bitbucket.org/deboghosh/hypar

  It can be cloned using git as follows:
  + git clone git@bitbucket.org:deboghosh/hypar.git (if you have a Bitbucket account)
  + git clone https://bitbucket.org/deboghosh/hypar.git (if you don't have a Bitbucket account)

  Bitbucket also allows downloading the package as a tarball, see 
  https://bitbucket.org/deboghosh/hypar/downloads.

  Documentation
  -------------
  To generate a local copy of this documentation, run "doxygen Doxyfile" in $(root_dir). The folder $(root_dir)/doc
  should contain the generated documentation in HTML format.

  Compiling
  ---------

  To compile HyPar, follow these steps in the root directory:
  
        autoreconf -i
        [CFLAGS="..."] ./configure [options]
        make
        make install

  CFLAGS should include all the compiler flags.

  \b Note: Default installation target is its own directory, and thus "make install" should not require
           administrative privileges. The binary will be placed in \a bin/ subdirectory.

  The configure options can include options such as BLAS/LAPACK location, MPI directory, etc. Type "./configure --help"
  to see a full list. The options specific to HyPar are:
  + --with-mpi-dir: Specify path where mpicc is installed.
  + --enable-omp: Enable OpenMP threads.
  + --enable-scalapack: Enable ScaLAPACK (this will make available a tridiagonal solver using ScaLAPACK).
  + --with-blas-dir: Specify path where BLAS is installed (relevant only if --enable-scalapack is specified).
  + --with-lapack-dir: Specify path where LAPACK is installed (relevant only if --enable-scalapack is specified).
  + --with-scalapack-dir: Specify path where ScaLAPACK is installed (relevant only if --enable-scalapack is specified).
  + --with-fortran-lib: Specify path where FORTRAN libraries are installed (for ScaLAPACK) (relevant only if --enable-scalapack 
    is specified).

  \b Compiling \b with \b PETSc:
  Install PETSc and make sure the environment variables \b PETSC_DIR and \b PETSC_ARCH are defined. Please PETSc's
  installation instructions for this. Once these environment variables are present, HyPar will use them to compile
  itself with PETSc functionalities.

  Notes
  -----
  + This package has been tested using the GNU and IBM C compilers. The configuration script is designed to look for these 
    compilers only.
  + Feel free to contact me about anything regarding this (doubts/difficulties/suggestions).
  + Feel free to use and modify the code in any way.

  Running
  -------
  + It's best to start with some examples. See the section on examples.
  + To run more cases, see the section in input files for a complete description of input files required.

  Testing and Baselines
  ---------------------
  + A set of baselines are available at: https://bitbucket.org/deboghosh/hypar_baselines
    These were obtained using the master branch of the main HyPar repo.
  + The script Extras/generate_baselines.sh can be used to clone this baselines repo and
    regenerate the solutions.
  + The scripts Extras/test_local.sh and Extras/test_repo.sh can be used to test a local or
    remote copy (git repo and branch) of HyPar against these baselines.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#ifdef with_petsc
#include <petscinterface.h>
#endif
#include <simulation.h>

static const char help[] = "HyPar - A finite-difference algorithm for solving hyperbolic-parabolic PDEs";

/*! \brief Initialize simulation objects
 *
 * Read in the number of simulations and allocate array of these objects.
*/
int InitializeSimulation( SimulationObject**  sim,    /*!< Array of simulation objects of type 
                                                           #SimulationObject, must be NULL. */
                          int*                nsims,  /*!< Number of simulation objects */
                          int                 rank,   /*!< MPI rank of this process */
                          int                 nproc   /*!< Number of MPI processes  */
                        )
{
  if (*sim != NULL) {
    fprintf(stderr,"Errror: sim is not NULL on rank %d.\n",rank);
    return 1;
  }

  /* default value */
  *nsims = 1;

  if (!rank) {

    FILE *in;
    in = fopen("simulation.inp","r");
    if (in) {
      int ferr;
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
        while (strcmp(word, "end")) {
  	      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
          if (!strcmp(word, "nsims")) {
            ferr = fscanf(in,"%d",nsims); if (ferr != 1) return(1);
          } else if (strcmp(word, "end")) {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless);
            printf("Warning: keyword %s in file \"simulation.inp\" with value %s not recognized or extraneous. Ignoring.\n",
                    word,useless);
          }
          if (ferr != 1) return(1);
        }
      } else {
   		  fprintf(stderr,"Error: Illegal format in file \"solver.inp\".\n");
        return 1;
      }
      fclose(in);
    }
    if (*nsims < 1) {
      fprintf(stderr,"Error in InitializeSimulation(): invalid value for nsims (%d)!\n", nsims);
      return 1;
    }
    printf("Number of simulation domains: %d\n", *nsims);
  }

#ifndef serial
  MPI_Bcast(nsims,1,MPI_INT,0,MPI_COMM_WORLD);
#endif

  int ns;
  SimulationObject* sim_array = (SimulationObject*) calloc( *nsims, sizeof(SimulationObject));
  for (ns = 0; ns < *nsims; ns++) {
    sim_array[ns].solver.my_idx = ns;
    sim_array[ns].solver.nsims = *nsims;
    sim_array[ns].mpi.rank = rank;
    sim_array[ns].mpi.nproc = nproc;
  }

  *sim = sim_array;

  if (!rank) {
    printf("Allocated simulation object(s).\n");
  }

  return 0;
}

/*!
 * \brief Main driver
 *
 * The main driver function that calls the initialization, solving, and cleaning up functions.
*/
int main(int argc,char **argv)
{
  SimulationObject  *sim;
  int               ierr = 0, d, n;
  struct timeval    main_start, solve_start;
  struct timeval    main_end  , solve_end  ;
#ifdef with_petsc
  PetscBool         use_petscts;
#endif
  int               use_petsc = 0;

  int rank, nproc;
  MPI_Comm world;
#ifdef serial
  rank  = 0;
  nproc = 1;
  printf("HyPar - Serial Version\n");
#else
  MPI_Init(&argc,&argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &world);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank );
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  if (!rank) printf("HyPar - Parallel (MPI) version with %d processes\n",nproc);
#endif

#ifdef with_petsc
  PetscInitialize(&argc,&argv,(char*)0,help);
  if (!rank) printf("Compiled with PETSc time integration.\n");
#endif

  gettimeofday(&main_start,NULL);

  /* Allocate simulation objects */
  int nsims;
  sim = NULL;
  ierr = InitializeSimulation(&sim, &nsims, rank, nproc);
  if (sim == NULL) {
    printf("Error: InitializeSimulation() failed to allocate simulation objects on rank %d\n",
           rank);
    return 1;
  }
  if (ierr) {
    printf("Error: InitializeSimulation() returned with status %d on process %d.\n",
            ierr,rank);
    return(ierr);
  }

#ifndef serial
  for (n = 0; n < nsims; n++) {
    MPI_Comm_dup(MPI_COMM_WORLD, &(sim[n].mpi.world));
  }
#endif
#ifdef with_petsc
  use_petscts = PETSC_FALSE; /* default value */
  ierr = PetscOptionsGetBool( PETSC_NULL,PETSC_NULL,
                              "-use-petscts",
                              &use_petscts,
                              PETSC_NULL); CHKERRQ(ierr);
  if (use_petscts == PETSC_TRUE) use_petsc = 1;
  for (n = 0; n < nsims; n++) {
    sim[n].solver.use_petscTS  = use_petscts;
  }
#endif

  /* Read Inputs */
  ierr = ReadInputs(sim, nsims, rank);
  if (ierr) {
    printf("Error: ReadInputs() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
  
  /* Initialize and allocate arrays */
  ierr = Initialize(sim, nsims);
  if (ierr) {
    printf("Error: Initialize() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
  
  /* read and set grid & initial solution */
  ierr = InitialSolution(sim, nsims);
  if (ierr) {
    printf("Error: InitialSolution() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
  
  /* Initialize domain boundaries */
  ierr = InitializeBoundaries(sim, nsims);
  if (ierr) {
    printf("Error: InitializeBoundaries() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
  
  /* Initialize immersed boundaries */
  ierr = InitializeImmersedBoundaries(sim, nsims);
  if (ierr) {
    printf("Error: InitializeImmersedBoundaries() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
  
  /* Initialize solvers */
  ierr = InitializeSolvers(sim, nsims);
  if (ierr) {
    printf("Error: InitializeSolvers() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
  
  /* Initialize physics */
  ierr = InitializePhysics(sim, nsims);
  if (ierr) {
    printf("Error: InitializePhysics() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
  
  /* Initializations complete */
  
  /* Run the solver */
#ifndef serial
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  gettimeofday(&solve_start,NULL);
#ifdef with_petsc
  if (use_petsc == 1) {
    /* Use PETSc time-integration */
    ierr = SolvePETSc(sim, nsims, rank, nproc);
    if (ierr) {
      printf("Error: SolvePETSc() returned with status %d on process %d.\n",ierr,rank);
      return(ierr);
    }
  } else {
    /* Use native time-integration */
    ierr = Solve(sim, nsims, rank, nproc);
    if (ierr) {
      printf("Error: Solve() returned with status %d on process %d.\n",ierr,rank);
      return(ierr);
    }
  }
#else 
  /* Use native time-integration */
  ierr = Solve(sim, nsims, rank, nproc);
  if (ierr) {
    printf("Error: Solve() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
#endif
  gettimeofday(&solve_end,NULL);
#ifndef serial
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  gettimeofday(&main_end,NULL);

  /* calculate solver and total runtimes */
  long long walltime;
  walltime = (  (main_end.tv_sec * 1000000   + main_end.tv_usec  ) 
              - (main_start.tv_sec * 1000000 + main_start.tv_usec));
  double main_runtime = (double) walltime / 1000000.0;
  ierr = MPIMax_double(&main_runtime,&main_runtime,1,&world); if(ierr) return(ierr);
  walltime = (  (solve_end.tv_sec * 1000000   + solve_end.tv_usec  ) 
              - (solve_start.tv_sec * 1000000 + solve_start.tv_usec));
  double solver_runtime = (double) walltime / 1000000.0;
  ierr = MPIMax_double(&solver_runtime,&solver_runtime,1,&world); if(ierr) return(ierr);

  /* Write errors and other data */
  SimWriteErrors(sim, nsims, rank, solver_runtime, main_runtime);

  /* Cleaning up */
  ierr = Cleanup(sim, nsims);
  if (ierr) {
    printf("Error: CleanUp() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
  free(sim);
  if (!rank) printf("Finished.\n");

#ifdef with_petsc
  PetscFinalize();
#endif

#ifndef serial
  MPI_Comm_free(&world);
  MPI_Finalize();
#endif
  return(0);
}
