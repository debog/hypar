/*! @file SparseGridsMiscFuncs.cpp
    @brief Miscellaneous functions needed for a sparse grids simulation
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <cfloat>
#include <limits>
#include <math.h>
#include <sparse_grids_simulation.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <io_cpp.h>
#include <std_vec_ops.h>

#ifdef with_python
#include <Python.h>
#endif

extern "C" void IncrementFilenameIndex(char*,int);

/*! Compute the number, coefficients, and dimensions of all the sparse
 * grids that go into the combination technique.
*/
int SparseGridsSimulation::ComputeSGDimsAndCoeffs()
{
  m_n_fg = log2(m_sim_fg->solver.dim_global[0]);
  int d = m_ndims;

  double c1 = 1;
  m_combination.clear();
  for (int s = 1; s <= d; s++) {

    double coeff = c1 * cfunc(d-1,s-1);

    std::vector<GridDimensions> dims(0);
    GetCTGridSizes((m_n_fg+(m_ndims-1)*(m_imin-1))+d-s, dims);
    if (dims.size() == 0) {
      fprintf(stderr, "Error in SparseGridsSimulation::ComputeSGDimsAndCoeffs()\n");
      fprintf(stderr, "  GetCTGridSize() returned empty vector!\n");
      m_combination.clear();
      return 1;
    }
    for (int i=0; i<dims.size(); i++) {
      m_combination.push_back(std::pair<double,GridDimensions>(coeff,dims[i]));
    }

    c1 *= -1;

  }

  return 0;
}

/*! Compute all grids such that the base2 logs of the size in each dimension
 * adds up to the input argument */
void SparseGridsSimulation::GetCTGridSizes(const int a_N, /*!< Desired sum of the base2 logs of grid sizes */
                                          std::vector<GridDimensions>&  a_dims  /*!< Vector of grid sizes */
                                         )
{
  a_dims.clear();

  int index[m_ndims], ubounds[m_ndims], lbounds[m_ndims];
  _ArraySetValue_(index,  m_ndims, m_imin);
  _ArraySetValue_(ubounds, m_ndims, a_N);
  _ArraySetValue_(lbounds, m_ndims, m_imin);

  int done = 0;
  while (!done) {
    int sum; _ArraySum1D_(index, sum, m_ndims);
    if (sum == a_N) {
      std::vector<int> tmp(m_ndims);
      for (int d=0; d<m_ndims; d++) {
        raiseto_int( tmp[d], 2,index[d] );
      }
      a_dims.push_back(tmp);
    }
    _ArrayIncrementIndexWithLBound_(m_ndims,ubounds,lbounds,index,done);
  }

  return;
}

/*! Compute the load-balanced processor distribution for a given grid size
 * and the total number of MPI ranks available */
int SparseGridsSimulation::ComputeProcessorDistribution(ProcDistribution&     a_iprocs, /*!< Processor distribution to compute */
                                                        const GridDimensions& a_dim     /*!< Grid dimensions */)
{
  a_iprocs.resize(m_ndims);
  /* get the normal vector for the grid dimensions */
  std::vector<double> dvec;
  StdVecOps::createNormalVector(dvec, a_dim);

  /* calculate the maximum number of MPI ranks in each dimension */
  int max_procs[m_ndims], min_procs[m_ndims];
  for (int i=0; i<m_ndims; i++) {
    max_procs[i] = a_dim[i]/_MIN_GRID_PTS_PER_PROC_;
    min_procs[i] = 1;
  }
  int max_nproc; _ArrayProduct1D_(max_procs, m_ndims, max_nproc);
  if (max_nproc < m_nproc) {
    fprintf(stderr, "Error in SparseGridsSimulation::ComputeProcessorDistribution() - rank %d\n", m_rank);
    fprintf(stderr, "  Number of MPI ranks greater than the maximum number of MPI ranks that can be used.\n");
    fprintf(stderr, "  Please re-run with %d MPI ranks.\n", max_nproc);
    return 1;
  }

  /* find all the processor distributions that are okay, i.e., their product
   * is the total number of MPI ranks #SparseGridsSimulation::m_nproc */
  std::vector<ProcDistribution> iproc_candidates(0);
  int iproc[m_ndims], ubound[m_ndims], done = 0;
  for (int d=0; d<m_ndims; d++) ubound[d] = max_procs[d]+1;
  _ArraySetValue_(iproc, m_ndims, 1);
  while (!done) {
    int prod; _ArrayProduct1D_(iproc, m_ndims, prod);
    if (prod == m_nproc) {
      ProcDistribution iproc_vec(m_ndims);
      for (int d = 0; d < m_ndims; d++) iproc_vec[d] = iproc[d];
      iproc_candidates.push_back(iproc_vec);
    }
    _ArrayIncrementIndexWithLBound_(m_ndims,ubound,min_procs,iproc,done);
  }
  if (iproc_candidates.size() == 0) {
    fprintf(stderr, "Error in SparseGridsSimulation::ComputeProcessorDistribution() - rank %d\n", m_rank);
    fprintf(stderr, "  Unable to fine candidate iprocs!\n");
    return 1;
  }

  /* find the candidate that is closest to the normalized dim vector */
  double min_norm = DBL_MAX;
  for (int i = 0; i<iproc_candidates.size(); i++) {
    std::vector<double> pvec;
    StdVecOps::createNormalVector(pvec, iproc_candidates[i]);
    double norm = StdVecOps::compute2Norm(pvec, dvec);
    if (norm < min_norm) {
      min_norm = norm;
      a_iprocs = iproc_candidates[i];
    }
  }

  return 0;
}

/*! Cleanup (deallocation) function for a  barebones simulation object that
    has been allocated using #SparseGridsSimulation::InitializeBarebones().
*/
int SparseGridsSimulation::CleanupBarebones( SimulationObject *sim /*!< simulation objects of type #SimulationObject*/)
{
  if (sim->is_barebones == 0) {
    fprintf(stderr, "Error in SparseGridsSimulation::CleanupBarebones()\n");
    fprintf(stderr, "  Simulation object is not a barebones one.\n");
    return 1;
  }

  HyPar* solver = &(sim->solver);
  MPIVariables* mpi = &(sim->mpi);

  /* Free the communicators created */
  IERR MPIFreeCommunicators(solver->ndims,mpi); CHECKERR(ierr);

  /* These variables are allocated in Initialize.c */
  free(solver->dim_global);
  free(solver->dim_global_ex);
  free(solver->dim_local);
  free(solver->index);
  free(solver->isPeriodic);
  free(solver->u);
  free(solver->x);
  free(solver->dxinv);
  free(mpi->iproc);
  free(mpi->ip);
  free(mpi->is);
  free(mpi->ie);
  free(mpi->bcperiodic);
  free(mpi->sendbuf);
  free(mpi->recvbuf);
  free(solver->VolumeIntegral);
  free(solver->VolumeIntegralInitial);
  free(solver->TotalBoundaryIntegral);
  free(solver->ConservationError);
  free(solver->stride_with_ghosts);
  free(solver->stride_without_ghosts);
  if (solver->filename_index) free(solver->filename_index);

  return(0);
}

/*! A barebones initialization function for a simulation object. It will
    allocate data holders for only the stuff needed for storing a solution;
    this object *cannot* be used for an actual simulation.
    This function does the following:
    + initializes the values for MPI variables
    + allocates memory for arrays to store full grid solution
*/
int SparseGridsSimulation::InitializeBarebones( SimulationObject *simobj /*!< simulation objects of type #SimulationObject*/)
{
  /* Set this to "true" since this is a barebones initialization */
  simobj->is_barebones = 1;

  HyPar* solver = &(simobj->solver);
  MPIVariables* mpi = &(simobj->mpi);

  /* allocations */
  mpi->ip             = (int*) calloc (m_ndims,sizeof(int));
  mpi->is             = (int*) calloc (m_ndims,sizeof(int));
  mpi->ie             = (int*) calloc (m_ndims,sizeof(int));
  mpi->bcperiodic     = (int*) calloc (m_ndims,sizeof(int));
  solver->dim_local   = (int*) calloc (m_ndims,sizeof(int));
  solver->isPeriodic  = (int*) calloc (m_ndims,sizeof(int));

#ifndef serial
  /* Domain partitioning */
  int total_proc = 1;
  for (int i=0; i<m_ndims; i++) total_proc *= mpi->iproc[i];

  /* calculate ndims-D rank of each process (ip[]) from rank in MPI_COMM_WORLD */
  IERR MPIRanknD( m_ndims,
                  m_rank,
                  mpi->iproc,
                  mpi->ip ); CHECKERR(ierr);

  /* calculate local domain sizes along each dimension */
  for (int i=0; i<m_ndims; i++) {
    solver->dim_local[i] = MPIPartition1D( solver->dim_global[i],
                                           mpi->iproc[i],
                                           mpi->ip[i] );
  }

  /* calculate local domain limits in terms of global domain */
  IERR MPILocalDomainLimits(  m_ndims,
                              m_rank,
                              &(simobj->mpi),
                              solver->dim_global,
                              mpi->is,
                              mpi->ie  );
  CHECKERR(ierr);

  /* create sub-communicators for parallel computations
   * along grid lines in each dimension */
   IERR MPICreateCommunicators(m_ndims,&(simobj->mpi)); CHECKERR(ierr);

  /* initialize periodic BC flags to zero */
  for (int i=0; i<solver->ndims; i++) mpi->bcperiodic[i] = 0;

  /* create communication groups */
  IERR MPICreateIOGroups(&(simobj->mpi)); CHECKERR(ierr);

#else

  for (int i=0; i<m_ndims; i++) {
    mpi->ip[i]            = 0;
    solver->dim_local[i]  = solver->dim_global[i];
    mpi->iproc[i]         = 1;
    mpi->is[i]            = 0;
    mpi->ie[i]            = solver->dim_local[i];
    mpi->bcperiodic[i]    = 0;
  }

#endif

  solver->npoints_global
    = solver->npoints_local
    = solver->npoints_local_wghosts
    = 1;
  for (int i=0; i<m_ndims; i++) {
    solver->npoints_global *= solver->dim_global[i];
    solver->npoints_local *= solver->dim_local [i];
    solver->npoints_local_wghosts *= (solver->dim_local[i]+2*solver->ghosts);
  }

  /* Allocations */
  if (!m_rank) printf("Allocating data arrays for full grid.\n");
  solver->index = (int*) calloc (m_ndims,sizeof(int));
  solver->stride_with_ghosts = (int*) calloc (solver->ndims,sizeof(int));
  solver->stride_without_ghosts = (int*) calloc (solver->ndims,sizeof(int));
  int accu1 = 1, accu2 = 1;
  for (int i=0; i<solver->ndims; i++) {
    solver->stride_with_ghosts[i]    = accu1;
    solver->stride_without_ghosts[i] = accu2;
    accu1 *= (solver->dim_local[i]+2*solver->ghosts);
    accu2 *=  solver->dim_local[i];
  }


  int size;
  /* state variables */
  size = 1;
  for (int i=0; i<m_ndims; i++) size *= (solver->dim_local[i]+2*solver->ghosts);
  solver->u = (double*) calloc (solver->nvars*size,sizeof(double));

  /* grid */
  size = 0;
  for (int i=0; i<m_ndims; i++) size += (solver->dim_local[i]+2*solver->ghosts);
  solver->x     = (double*) calloc (size,sizeof(double));
  solver->dxinv = (double*) calloc (size,sizeof(double));

  /* allocate MPI send/receive buffer arrays */
  int bufdim[solver->ndims], maxbuf = 0;
  for (int d = 0; d < solver->ndims; d++) {
    bufdim[d] = 1;
    for (int i = 0; i < solver->ndims; i++) {
      if (i == d) bufdim[d] *= solver->ghosts;
      else        bufdim[d] *= solver->dim_local[i];
    }
    if (bufdim[d] > maxbuf) maxbuf = bufdim[d];
  }
  maxbuf *= solver->nvars;
  mpi->maxbuf  = maxbuf;
  mpi->sendbuf = (double*) calloc (2*solver->ndims*maxbuf,sizeof(double));
  mpi->recvbuf = (double*) calloc (2*solver->ndims*maxbuf,sizeof(double));

  solver->VolumeIntegral        = (double*) calloc (solver->nvars,sizeof(double));
  solver->VolumeIntegralInitial = (double*) calloc (solver->nvars,sizeof(double));
  solver->TotalBoundaryIntegral = (double*) calloc (solver->nvars,sizeof(double));
  solver->ConservationError     = (double*) calloc (solver->nvars,sizeof(double));

  for (int i=0; i<solver->nvars; i++) solver->ConservationError[i] = -1;

  return(0);
}

/*! This function initializes all some function pointers needed by a
 *  barebones simulation object.
*/
int SparseGridsSimulation::InitializeSolversBarebones( SimulationObject *sim /*!< simulation objects of type #SimulationObject*/)
{
  if (sim->is_barebones != 1) {
    fprintf(stderr, "Error in SparseGridsSimulation::InitializeSolversBarebones() - \n");
    fprintf(stderr, "  simulation object is not barebones type.\n");
  }

  int ns;

  HyPar        *solver   = &(sim->solver);
  MPIVariables *mpi      = &(sim->mpi);

  solver->ParabolicFunction         = NULL;
  solver->SecondDerivativePar       = NULL;
  solver->FirstDerivativePar        = NULL;
  solver->InterpolateInterfacesPar  = NULL;

  solver->interp                = NULL;
  solver->compact               = NULL;
  solver->lusolver              = NULL;
  solver->SetInterpLimiterVar   = NULL;
  solver->flag_nonlinearinterp  = 0;
  solver->time_integrator = NULL;
  solver->msti = NULL;

  /* Solution output function */
  solver->WriteOutput    = NULL; /* default - no output */
  solver->filename_index = NULL;
  if (!strcmp(solver->output_mode,"serial")) {
    solver->index_length = 5;
    solver->filename_index = (char*) calloc (solver->index_length+1,sizeof(char));
    int i; for (i=0; i<solver->index_length; i++) solver->filename_index[i] = '0';
    solver->filename_index[solver->index_length] = (char) 0;
    if (!strcmp(solver->op_file_format,"text")) {
      solver->WriteOutput = WriteText;
      strcpy(solver->solnfilename_extn,".dat");
    } else if (!strcmp(solver->op_file_format,"tecplot2d")) {
      solver->WriteOutput = WriteTecplot2D;
      strcpy(solver->solnfilename_extn,".dat");
    } else if (!strcmp(solver->op_file_format,"tecplot3d")) {
      solver->WriteOutput = WriteTecplot3D;
      strcpy(solver->solnfilename_extn,".dat");
    } else if ((!strcmp(solver->op_file_format,"binary")) || (!strcmp(solver->op_file_format,"bin"))) {
      solver->WriteOutput = WriteBinary;
      strcpy(solver->solnfilename_extn,".bin");
    } else if (!strcmp(solver->op_file_format,"none")) {
      solver->WriteOutput = NULL;
    } else {
      fprintf(stderr,"Error (domain %d): %s is not a supported file format.\n",
              ns, solver->op_file_format);
      return(1);
    }
    if ((!strcmp(solver->op_overwrite,"no")) && solver->restart_iter) {
      /* if it's a restart run, fast-forward the filename */
      int t;
      for (t=0; t<solver->restart_iter; t++)
        if ((t+1)%solver->file_op_iter == 0) IncrementFilenameIndex(solver->filename_index,solver->index_length);
    }
  } else if (!strcmp(solver->output_mode,"parallel")) {
    if (!strcmp(solver->op_file_format,"none")) solver->WriteOutput = NULL;
    else {
      /* only binary file writing supported in parallel mode */
      /* use post-processing scripts to convert              */
      solver->WriteOutput = WriteBinary;
      strcpy(solver->solnfilename_extn,".bin");
    }
  } else {

    fprintf(stderr,"Error (domain %d): %s is not a supported output mode.\n",
            ns, solver->output_mode);
    fprintf(stderr,"Should be \"serial\" or \"parallel\".    \n");
    return(1);

  }

  /* Solution plotting function */
  strcpy(solver->plotfilename_extn,".png");
#ifdef with_python
  solver->py_plt_func = NULL;
  solver->py_plt_func_args = NULL;
  {
    char python_plotting_fname[_MAX_STRING_SIZE_] = "plotSolution";
    PyObject* py_plot_name = PyUnicode_DecodeFSDefault(python_plotting_fname);
    PyObject* py_plot_module = PyImport_Import(py_plot_name);
    Py_DECREF(py_plot_name);
    if (py_plot_module) {
      solver->py_plt_func = PyObject_GetAttrString(py_plot_module, "plotSolution");
      if (!solver->py_plt_func) {
        if (!mpi->rank) {
          printf("Unable to load plotSolution function from Python module.\n");
        }
      } else {
        if (!mpi->rank) {
          printf("Loaded Python module for plotting.\n");
          printf("Loaded plotSolution function from Python module.\n");
        }
      }
    } else {
      if (!mpi->rank) {
        printf("Unable to load Python module for plotting.\n");
      }
    }
  }
#endif

  return 0;
}
