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
  + Can use sparse grids for faster computations on high-dimensional problems
  + Written entirely in C/C++ and uses the MPICH library. It also uses OpenMP threads 
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
        [CFLAGS="..."] [CXXFLAGS="..."] ./configure [options]
        make
        make install

  CFLAGS and CXXFLAGS should include all the compiler flags.

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
  + This package has been tested using the GNU C and C++ compilers. The configuration script is designed to look for these 
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
#include <sys/time.h>
#include <string>

#ifdef with_petsc
#include <petscinterface.h>
#endif
#include <mpivars_cpp.h>
#include <simulation_library.h>

static const char help[] = "HyPar - A finite-difference algorithm for solving hyperbolic-parabolic PDEs";

/*!
 * \brief Main driver
 *
 * The main driver function that calls the initialization, solving, and cleaning up functions.
*/
int main(int argc, char **argv)
{
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

  int sim_type = -1;
  Simulation *sim = NULL;

  if (!rank) {

    std::string ensemble_sim_fname(_ENSEMBLE_SIM_INP_FNAME_);
    std::string sparsegrids_sim_fname(_SPARSEGRIDS_SIM_INP_FNAME_);

    FILE *f_ensemble_sim = fopen(ensemble_sim_fname.c_str(), "r");
    FILE *f_sparsegrids_sim = fopen(sparsegrids_sim_fname.c_str(), "r");

    if (f_ensemble_sim && f_sparsegrids_sim) {

      fprintf(stderr,"Error: Cannot have both %s and %s input files.\n",
              _ENSEMBLE_SIM_INP_FNAME_, _SPARSEGRIDS_SIM_INP_FNAME_);
      fprintf(stderr, "Remove one or both of them depending on the kind of simulation you want to run.\n");
      fclose(f_ensemble_sim);
      fclose(f_sparsegrids_sim);

    } else if (f_ensemble_sim) {

      sim_type = _SIM_TYPE_ENSEMBLE_;
      fclose(f_ensemble_sim);

    } else if (f_sparsegrids_sim) {

      sim_type = _SIM_TYPE_SPARSE_GRIDS_;
      fclose(f_sparsegrids_sim);

    } else {

      sim_type = _SIM_TYPE_SINGLE_;

    }

  }

#ifndef serial
  MPI_Bcast(&sim_type, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

  if (sim_type == _SIM_TYPE_SINGLE_) {
    sim = new SingleSimulation;
  } else if (sim_type == _SIM_TYPE_ENSEMBLE_) {
    if (!rank) printf("-- Ensemble Simulation --\n");
    sim = new EnsembleSimulation;
  } else if (sim_type == _SIM_TYPE_SPARSE_GRIDS_) {
    if (!rank) printf("-- Sparse Grids Simulation --\n");
    sim = new SparseGridsSimulation;
  } else {
    fprintf(stderr, "ERROR: invalid sim_type (%d) on rank %d.\n",
            sim_type, rank);
  }

  if (sim == NULL) {
    fprintf(stderr, "ERROR: unable to create sim on rank %d.\n",
            rank );
    return 1;
  }

  /* Allocate simulation objects */
  ierr = sim->define(rank, nproc);
  if (!sim->isDefined()) {
    printf("Error: Simulation::define() failed on rank %d\n",
           rank);
    return 1;
  }
  if (ierr) {
    printf("Error: Simulation::define() returned with status %d on process %d.\n",
            ierr, rank);
    return(ierr);
  }

#ifndef serial
  ierr = sim->mpiCommDup();
#endif

#ifdef with_petsc
  use_petscts = PETSC_FALSE; /* default value */
  ierr = PetscOptionsGetBool( PETSC_NULL,PETSC_NULL,
                              "-use-petscts",
                              &use_petscts,
                              PETSC_NULL); CHKERRQ(ierr);
  if (use_petscts == PETSC_TRUE) use_petsc = 1;
  sim->usePetscTS(use_petscts);
#endif

  /* Read Inputs */
  ierr = sim->ReadInputs();
  if (ierr) {
    printf("Error: Simulation::ReadInputs() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
  
  /* Initialize and allocate arrays */
  ierr = sim->Initialize();
  if (ierr) {
    printf("Error: Simulation::Initialize() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
  
  /* read and set grid & initial solution */
  ierr = sim->InitialSolution();
  if (ierr) {
    printf("Error: Simulation::InitialSolution() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
  
  /* Initialize domain boundaries */
  ierr = sim->InitializeBoundaries();
  if (ierr) {
    printf("Error: Simulation::InitializeBoundaries() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
  
  /* Initialize immersed boundaries */
  ierr = sim->InitializeImmersedBoundaries();
  if (ierr) {
    printf("Error: Simulation::InitializeImmersedBoundaries() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
  
  /* Initialize solvers */
  ierr = sim->InitializeSolvers();
  if (ierr) {
    printf("Error: Simulation::InitializeSolvers() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }
  
  /* Initialize physics */
  ierr = sim->InitializePhysics();
  if (ierr) {
    printf("Error: Simulation::InitializePhysics() returned with status %d on process %d.\n",ierr,rank);
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
    ierr = sim->SolvePETSc();
    if (ierr) {
      printf("Error: Simulation::SolvePETSc() returned with status %d on process %d.\n",ierr,rank);
      return(ierr);
    }
  } else {
    /* Use native time-integration */
    ierr = sim->Solve();
    if (ierr) {
      printf("Error: Simulation::Solve() returned with status %d on process %d.\n",ierr,rank);
      return(ierr);
    }
  }
#else 
  /* Use native time-integration */
  ierr = sim->Solve();
  if (ierr) {
    printf("Error: Simulation::Solve() returned with status %d on process %d.\n",ierr,rank);
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
  sim->WriteErrors(solver_runtime, main_runtime);

  /* Cleaning up */
  delete sim;
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
