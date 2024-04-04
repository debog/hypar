/*! @file SparseGridsCalculateError.cpp
    @brief Calculate errors if an exact solution is provided
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <arrayfunctions.h>
#include <mathfunctions_cpp.h>
#include <std_vec_ops.h>
#include <sparse_grids_simulation.h>

extern "C" int ExactSolution(void*,void*,double*,char*,int*);

/*! Calculates the error in the solution if the exact solution is
    available. The exact solution should be provided in the file
    "exact.inp" in the same format as the initial solution.
*/
void SparseGridsSimulation::CalculateError()
{
  int exact_flag;
  double *uex = NULL;

  /* read in full grid exact solution, if available */
  {
    HyPar* solver = &(m_sim_fg->solver);
    MPIVariables* mpi = &(m_sim_fg->mpi);
    long size = solver->nvars * solver->npoints_local_wghosts;
    uex = (double*) calloc (size, sizeof(double));

    char fname_root[_MAX_STRING_SIZE_] = "exact";
    IERR ExactSolution( solver,
                        mpi,
                        uex,
                        fname_root,
                        &exact_flag ); CHECKERR(ierr);
  }

  for (int n=0; n<m_nsims_sg; n++) {
//    TimeError(  &(m_sims_sg[n].solver),
//                &(m_sims_sg[n].mpi),
//                uex );
  }

  if (!exact_flag) {

    /* No exact solution available */
    for (int n=0; n<m_nsims_sg; n++) {
      m_sims_sg[n].solver.error[0]
        = m_sims_sg[n].solver.error[1]
        = m_sims_sg[n].solver.error[2]
        = -1;
    }
    m_sim_fg->solver.error[0]
      = m_sim_fg->solver.error[1]
      = m_sim_fg->solver.error[2]
      = -1;

  } else {

    std::vector<int> periodic_arr(m_ndims);
    for (int i=0; i<m_ndims; i++) {
      periodic_arr[i] = (m_is_periodic[i] ? 1 : 0);
    }

    double *uex2 = NULL;

    if (m_print_sg_errors == 1) {
      long size =   m_sim_fg->solver.nvars
                  * m_sim_fg->solver.npoints_local_wghosts;
      uex2 = (double*) calloc(size, sizeof(double));
      _ArrayCopy1D_(uex, uex2, size);
    }

    /* calculate error for full grid */
    computeError( *m_sim_fg, uex );

    /* calculate error for sparse grids */
    if (m_print_sg_errors == 1) {
      for (int n = 0; n < m_nsims_sg; n++) {

        GridDimensions dim_fg(m_ndims,0);
        StdVecOps::copyFrom(dim_fg, m_sim_fg->solver.dim_global, m_ndims);

        GridDimensions dim_sg(m_ndims,0);
        StdVecOps::copyFrom(dim_sg, m_sims_sg[n].solver.dim_global, m_ndims);

        /* assemble the global exact solution on full grid */
        double *uex_global_fg = NULL;
        if (!m_rank) {
          allocateDataArrays( dim_fg,
                              m_sim_fg->solver.nvars,
                              &uex_global_fg,
                              m_sim_fg->solver.ghosts);
        }
        MPIGatherArraynDwGhosts( m_ndims,
                                 (void*) &(m_sim_fg->mpi),
                                 uex_global_fg,
                                 uex2,
                                 m_sim_fg->solver.dim_global,
                                 m_sim_fg->solver.dim_local,
                                 m_sim_fg->solver.ghosts,
                                 m_sim_fg->solver.nvars );

        /* interpolate to sparse grid -
         * this will delete the full grid array*/
        double *uex_global_sg = NULL;
        if (!m_rank) {
          int ierr = ::InterpolateGlobalnDVar(  dim_sg.data(),
                                                &uex_global_sg,
                                                dim_fg.data(),
                                                uex_global_fg,
                                                m_sims_sg[n].solver.nvars,
                                                m_sims_sg[n].solver.ghosts,
                                                m_ndims,
                                                periodic_arr.data() );
          if (ierr) {
            fprintf(stderr,"InterpolateGlobalnDVar() returned with error!\n");
            exit(1);
          }
        }

        /* allocate local exact solution on this sparse grid */
        long size = m_sims_sg[n].solver.nvars
                    * m_sims_sg[n].solver.npoints_local_wghosts;
        double* uex_sg = (double*) calloc(size, sizeof(double));

        /* partition the global exact solution to local on this sparse grid */
        MPIPartitionArraynDwGhosts( m_ndims,
                                    (void*) &(m_sims_sg[n].mpi),
                                    (m_rank ? NULL : uex_global_sg),
                                    uex_sg,
                                    m_sims_sg[n].solver.dim_global,
                                    m_sims_sg[n].solver.dim_local,
                                    m_sims_sg[n].solver.ghosts,
                                    m_sims_sg[n].solver.nvars );

        /* delete the global exact solution array */
        if (!m_rank) free(uex_global_sg);

        /* compute the error */
        computeError( m_sims_sg[n], uex_sg);

        /* delete the exact solution array */
        free(uex_sg);
      }
      free(uex2);
    }

  }

  free(uex);
  return;
}

/*! compute errors for a particular simulation object */
void SparseGridsSimulation::computeError( SimulationObject& a_sim,  /*!< Simulation object */
                                          double* const     a_uex   /*!< Exact solution */ )
{
  static const double tolerance = 1e-15;

  HyPar* solver = &(a_sim.solver);
  MPIVariables* mpi = &(a_sim.mpi);

  if (a_uex == NULL) {
    fprintf(stderr, "Error in SparseGrids::computeError() -\n");
    fprintf(stderr, "  exact solution pointer in NULL.\n");
    exit(1);
  }

  /* calculate solution norms (for relative error) */
  double sum, global_sum;
  double solution_norm[3] = {0.0,0.0,0.0};
  /* L1 */
  sum = ArraySumAbsnD   (solver->nvars,solver->ndims,solver->dim_local,
                         solver->ghosts,solver->index,a_uex);
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
  solution_norm[0] = global_sum/((double)solver->npoints_global);
  /* L2 */
  sum = ArraySumSquarenD(solver->nvars,solver->ndims,solver->dim_local,
                         solver->ghosts,solver->index,a_uex);
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
  solution_norm[1] = sqrt(global_sum/((double)solver->npoints_global));
  /* Linf */
  sum = ArrayMaxnD      (solver->nvars,solver->ndims,solver->dim_local,
                         solver->ghosts,solver->index,a_uex);
  global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->world);
  solution_norm[2] = global_sum;

  /* compute error = difference between exact and numerical solution */
  long size = solver->nvars*solver->npoints_local_wghosts;
  _ArrayAXPY_(solver->u,-1.0,a_uex,size);

  /* calculate L1 norm of error */
  sum = ArraySumAbsnD   (solver->nvars,solver->ndims,solver->dim_local,
                         solver->ghosts,solver->index,a_uex);
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
  solver->error[0] = global_sum/((double)solver->npoints_global);

  /* calculate L2 norm of error */
  sum = ArraySumSquarenD(solver->nvars,solver->ndims,solver->dim_local,
                         solver->ghosts,solver->index,a_uex);
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
  solver->error[1] = sqrt(global_sum/((double)solver->npoints_global));

  /* calculate Linf norm of error */
  sum = ArrayMaxnD      (solver->nvars,solver->ndims,solver->dim_local,
                         solver->ghosts,solver->index,a_uex);
  global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->world);
  solver->error[2] = global_sum;

  /*
    decide whether to normalize and report relative errors,
    or report absolute errors.
  */
  if (    (solution_norm[0] > tolerance)
      &&  (solution_norm[1] > tolerance)
      &&  (solution_norm[2] > tolerance) ) {
    solver->error[0] /= solution_norm[0];
    solver->error[1] /= solution_norm[1];
    solver->error[2] /= solution_norm[2];
  }

  return;
}

