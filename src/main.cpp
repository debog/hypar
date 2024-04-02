/*! @file main.c
 *  @brief Main driver.
 * The main driver function that calls the initialization, solving, and cleaning up functions.
 *  @author Debojyoti Ghosh
*/

/*! @mainpage

  @author Debojyoti Ghosh [\b Email: (first name) (dot) (last name) (at) gmail (dot) com, \b Website: http://debog.github.io/]
  @author John Loffeld
  @author Youngdae Kim

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
  + Can use <B>sparse grids</B> for faster computations on high-dimensional problems
  + Written entirely in C/C++ and uses the MPICH library. It can use OpenMP threads and CUDA
    on NVidia GPUs, but these are works-in-progress.
  + Can be <B>compiled with PETSc</B> (https://petsc.org/release/), if available, where
    it can use PETSc's time integration module TS (https://petsc.org/release/src/ts/).
  + For 3-dimensional simulations, the <B>immersed boundaries</B> can be used to
    solve over non-Cartesian geometries.
  + Can be <B>compiled with libROM</B> (https://www.librom.net/), if available, where
    it can be used with the reduced order modeling tools implemented in libROM.
  + Can be <B>compiled with Python</B> (https://www.python.org/), if available, to enable
    features like in-situ plotting of solution.

  HyPar has been developed to be scalable, and apart from the usual functionalities to
  solve a system of PDEs on distributed architectures, it provides scalable file I/O
  functions. It has been tested on several platforms, including DOE Leadership-class
  supercomputers, with up to ~0.5 million MPI ranks.

  Download
  --------
  The code is available at:
  + https://github.com/debog/hypar
  + https://gitlab.com/debojyoti.ghosh/hypar

  Follow the instructions on these pages to clone using git. For example,
  + git clone git@gitlab.com:debojyoti.ghosh/hypar.git
  + git clone https://github.com/debog/hypar.git

  Both these platforms allow downloading the package as a zip/tarball, for example
  + https://github.com/debog/hypar/archive/refs/heads/master.zip
  + https://gitlab.com/debojyoti.ghosh/hypar/-/archive/master/hypar-master.zip

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
  + \--enable-serial: Compile a serial version without MPI.
  + \--with-mpi-dir: Specify path where mpicc is installed, if not in standard path.
  + \--enable-omp: Enable OpenMP threads.
  + \--enable-cuda: Enable CUDA if NVidia GPU present.
  + \--enable-python: Enable Python interface.
  + \--with-cuda-dir: Specify path where CUDA is installed, if not in standard path.
  + \--enable-scalapack: Enable ScaLAPACK (this will make available a tridiagonal solver using ScaLAPACK).
  + \--enable-fftw: Enable FFTW (this will make available features that use the FFTW library; *needs MPI*).
  + \--with-blas-dir: Specify path where BLAS is installed (relevant only if \--enable-scalapack is specified).
  + \--with-lapack-dir: Specify path where LAPACK is installed (relevant only if \--enable-scalapack is specified).
  + \--with-scalapack-dir: Specify path where ScaLAPACK is installed (relevant only if \--enable-scalapack is specified).
  + \--with-fftw-dir: Specify path where FFTW is installed (relevant only if \--enable-fftw is specified). Note that the FFTW library should be compiled with MPI.
  + \--with-fortran-lib: Specify path where FORTRAN libraries are installed (for ScaLAPACK) (relevant only if \--enable-scalapack
    is specified).

  \b Notes:

  + Limited parts of the code have been implemented on CUDA. See \b Examples for currently available
    CUDA-enabled simulations.
  + Interfacing with Python and the features it enables are a work-in-progress. It is not very robust; for example
    it may result in inexplicable segmentation faults depending on the Python version and the MPI library. It should
    work okay when compiled with OpenMPI and when compiled without any MPI (serial). It is based on the approach
    in PythonFOAM (https://github.com/argonne-lcf/PythonFOAM). The following environment variables must be
    consistently defined:
    + \b PYTHON_LIB_PATH - location of Python libraries
    + \b PYTHON_BIN_PATH - location of Python binary
    + \b PYTHON_INCLUDE_PATH - location of Python headers
    + \b NUMPY_INCLUDE_PATH - location of numpy headers (for eg., /path/to/python3.8/dist-packages/numpy/core/include)
    + \b PYTHON_LIB_NAME - name of the Python library (for eg., lpython3.8)

  Compiling with other scientific computing libraries
  ---------------------------------------------------

  HyPar can be compiled with and use the functionalities of the following libraries. It is assumed that the user
  is familiar with them and has access to their documentation.

  \b Compiling \b with \b PETSc:
  Install PETSc and make sure the environment variables \b PETSC_DIR and \b PETSC_ARCH are defined. Please see PETSc's
  installation instructions for this. Once these environment variables are present, HyPar will use them to compile
  itself with PETSc functionalities.

  \b Compiling \b with \b libROM:
  Download and compile libROM and make sure the environment variable \b LIBROM_DIR is defined as the location of libROM
  source directory. Once this environment variable is present, HyPar will use it to compile itself with libROM functionalities.

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
  + A set of baselines are available at: https://gitlab.com/debojyoti.ghosh/hypar_baselines
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

#ifdef with_python
#include <Python.h>
#ifdef with_python_numpy
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#endif
#endif

#include <mpivars_cpp.h>
#include <simulation_library.h>

static const char help[] = "HyPar - A finite-difference algorithm for solving hyperbolic-parabolic PDEs";

#ifdef with_python
static void initializePython(int);
static void initializePythonPlotting(int);
#endif

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

#ifdef serial
  int world = 0;
  int rank  = 0;
  int nproc = 1;
  printf("HyPar - Serial Version\n");
#else
  MPI_Comm world;
  int rank, nproc;
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

#ifdef with_python
  initializePython(rank);
  initializePythonPlotting(rank);
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
  ierr = PetscOptionsGetBool( nullptr,nullptr,
                              "-use-petscts",
                              &use_petscts,
                              nullptr); CHKERRQ(ierr);
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

  /* Initialize physics data */
  ierr = sim->InitializePhysicsData();
  if (ierr) {
    printf("Error: Simulation::InitializePhysicsData() returned with status %d on process %d.\n",ierr,rank);
    return(ierr);
  }

  /* Wrap up initializations */
  ierr = sim->InitializationWrapup();
  if (ierr) {
    printf("Error: Simulation::InitializationWrapup() returned with status %d on process %d.\n",ierr,rank);
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

#ifdef with_python
  Py_Finalize();
#endif

#ifdef with_petsc
  PetscFinalize();
#endif

#ifndef serial
  MPI_Comm_free(&world);
  MPI_Finalize();
#endif

  return(0);
}

#ifdef with_python
/*! Initialize Python interpreter */
void initializePython(int a_rank /*!< MPI rank */)
{
  Py_Initialize();
  if (!a_rank) printf("Initialized Python.\n");
  PyRun_SimpleString("import os");
  PyRun_SimpleString("hypar_dir = os.environ.get('HYPAR_DIR')");
#ifndef serial
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  return;
}

/*! Initialize Python plotting stuff */
void initializePythonPlotting(int a_rank /*!< MPI rank */)
{
  /* load plotting tools and scripts */
  PyRun_SimpleString("hypar_dir_plt_py = hypar_dir + '/src/PlottingFunctions'");
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append(hypar_dir_plt_py)");
  if (!a_rank) {
    PyRun_SimpleString
      ("print('Added plotting script directory (%s) to Python path.' % hypar_dir_plt_py)");
  }
#ifndef serial
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  return;
}

#endif
