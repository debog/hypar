/*! @file SparseGridsSolve.cpp
    @brief Run a sparse grids simulation
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <common_cpp.h>
#include <io_cpp.h>
#include <timeintegration_cpp.h>
#include <sparse_grids_simulation.h>

/*! This function integrates the semi-discrete ODE (obtained from discretizing the
    PDE in space) using natively implemented time integration methods. It initializes
    the time integration object, iterates the simulation for the required number of
    time steps, and calculates the errors. After the specified number of iterations,
    it writes out some information to the screen and the solution to a file.
*/
int SparseGridsSimulation::Solve()
{
  int tic     = 0;

  /* write out iblank to file for visualization */
  if (m_write_sg_solutions == 1) {
    for (int ns = 0; ns < m_nsims_sg; ns++) {
      if (m_sims_sg[ns].solver.flag_ib) {

        char fname_root[_MAX_STRING_SIZE_] = "iblank_sg";
        if (m_nsims_sg > 1) {
          char index[_MAX_STRING_SIZE_];
          GetStringFromInteger(ns, index, (int)log10((m_nsims_sg)+1));
          strcat(fname_root, "_");
          strcat(fname_root, index);
        }

        WriteArray( m_sims_sg[ns].solver.ndims,
                    1,
                    m_sims_sg[ns].solver.dim_global,
                    m_sims_sg[ns].solver.dim_local,
                    m_sims_sg[ns].solver.ghosts,
                    m_sims_sg[ns].solver.x,
                    m_sims_sg[ns].solver.iblank,
                    &(m_sims_sg[ns].solver),
                    &(m_sims_sg[ns].mpi),
                    fname_root );
      }
    }
  }
  if (m_sim_fg->solver.flag_ib) {

    char fname_root[_MAX_STRING_SIZE_] = "iblank";
    WriteArray( m_sim_fg->solver.ndims,
                1,
                m_sim_fg->solver.dim_global,
                m_sim_fg->solver.dim_local,
                m_sim_fg->solver.ghosts,
                m_sim_fg->solver.x,
                m_sim_fg->solver.iblank,
                &(m_sim_fg->solver),
                &(m_sim_fg->mpi),
                fname_root );
  }

  /* Define and initialize the time-integration object */
  TimeIntegration TS;
  if (!m_rank) printf("Setting up time integration.\n");
  TimeInitialize((void*)m_sims_sg.data(), m_nsims_sg, m_rank, m_nproc, &TS);

  if (!m_rank) {
    printf( "Solving in time (from %d to %d iterations)\n",
            TS.restart_iter,TS.n_iter);
  }

  for (TS.iter = TS.restart_iter; TS.iter < TS.n_iter; TS.iter++) {

    /* Write initial solution to file if this is the first iteration */
    if (!TS.iter) OutputSolution(TS.waqt);

    /* Call pre-step function */
    TimePreStep  (&TS);

    /* Step in time */
    TimeStep     (&TS);

    /* Call post-step function */
    TimePostStep (&TS);

    /* Print information to screen */
    TimePrintStep(&TS);
    tic++;

    /* Write intermediate solution to file */
    if ((TS.iter+1)%m_sims_sg[0].solver.file_op_iter == 0) {
      OutputSolution(TS.waqt);
      tic = 0;
    }

  }

  /* write a final solution file, if last iteration did not write one */
  if (tic || (!TS.n_iter)) OutputSolution(TS.waqt);

  if (!m_rank) {
    printf("Completed time integration (Final time: %f).\n",TS.waqt);
    if (m_nsims_sg > 1) printf("\n");
  }

  /* calculate error if exact solution has been provided */
  CalculateError();

  TimeCleanup(&TS);

  return(0);
}

