/*! @file Solve.c
    @author Debojyoti Ghosh
    @brief  Solve the governing equations in time
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <common_cpp.h>
#include <io_cpp.h>
#include <timeintegration_cpp.h>
#include <mpivars_cpp.h>
#include <simulation_object.h>

#ifdef with_librom
#include <librom_interface.h>
#endif

#ifdef compute_rhs_operators
extern "C" int ComputeRHSOperators(void*,void*,double);
#endif
extern "C" int CalculateError(void*,void*); /*!< Calculate the error in the final solution */
int OutputSolution(void*,int,double);   /*!< Write solutions to file */
extern "C" void ResetFilenameIndex(char*, int); /*!< Reset filename index */
#ifdef with_librom
extern "C" int CalculateROMDiff(void*,void*); /*!< Calculate the diff of PDE and ROM solutions */
int OutputROMSolution(void*,int,double);   /*!< Write ROM solutions to file */
#endif

/*! This function integrates the semi-discrete ODE (obtained from discretizing the
    PDE in space) using natively implemented time integration methods. It initializes
    the time integration object, iterates the simulation for the required number of
    time steps, and calculates the errors. After the specified number of iterations,
    it writes out some information to the screen and the solution to a file.
*/
int Solve(  void  *a_s,     /*!< Array of simulation objects of type #SimulationObject */
            int   a_nsims,  /*!< number of simulation objects */
            int   a_rank,   /*!< MPI a_rank of this process */
            int   nproc   /*!< Number of MPI processes */
         )
{
  SimulationObject* sim = (SimulationObject*) a_s;

  /* make sure none of the simulation objects sent in the array
   * are "barebones" type */
  for (int ns = 0; ns < a_nsims; ns++) {
    if (sim[ns].is_barebones == 1) {
      fprintf(stderr, "Error in Solve(): simulation object %d on a_rank %d is barebones!\n",
              ns, a_rank );
      return 1;
    }
  }

  /* write out iblank to file for visualization */
  for (int ns = 0; ns < a_nsims; ns++) {
    if (sim[ns].solver.m_flag_ib) {

      char fname_root[_MAX_STRING_SIZE_] = "iblank";
      if (a_nsims > 1) {
        char index[_MAX_STRING_SIZE_];
        GetStringFromInteger(ns, index, (int)log10((a_nsims)+1));
        strcat(fname_root, "_");
        strcat(fname_root, index);
      }

      WriteArray( sim[ns].solver.m_ndims,
                  1,
                  sim[ns].solver.m_dim_global,
                  sim[ns].solver.m_dim_local,
                  sim[ns].solver.m_ghosts,
                  sim[ns].solver.m_x,
                  sim[ns].solver.m_iblank,
                  &(sim[ns].solver),
                  &(sim[ns].mpi),
                  fname_root );
    }
  }

#ifdef with_librom
  if (!a_rank) printf("Setting up libROM interface.\n");
  libROMInterface rom_interface( sim, a_nsims, a_rank, nproc, sim[0].solver.m_dt );
  const std::string& rom_mode( rom_interface.mode() );
  std::vector<double> op_times_arr(0);
#endif

#ifdef with_librom
  if ((rom_mode == _ROM_MODE_TRAIN_) || (rom_mode == _ROM_MODE_NONE_)) {
#endif
    /* Define and initialize the time-integration object */
    TimeIntegration TS;
    if (!a_rank) printf("Setting up time integration.\n");
    TimeInitialize(sim, a_nsims, a_rank, nproc, &TS);
    double ti_runtime = 0.0;

    if (!a_rank) printf("Solving in time (from %d to %d iterations)\n",TS.m_restart_iter,TS.m_n_iter);
    for (TS.m_iter = TS.m_restart_iter; TS.m_iter < TS.m_n_iter; TS.m_iter++) {

      /* Write initial solution to file if this is the first iteration */
      if (!TS.m_iter) {
        for (int ns = 0; ns < a_nsims; ns++) {
          if (sim[ns].solver.PhysicsOutput) {
            sim[ns].solver.PhysicsOutput( &(sim[ns].solver),
                                          &(sim[ns].mpi),
                                          TS.m_waqt );
          }
        }
        OutputSolution(sim, a_nsims, TS.m_waqt);
#ifdef with_librom
        op_times_arr.push_back(TS.m_waqt);
#endif
      }

#ifdef with_librom
      if ((rom_mode == _ROM_MODE_TRAIN_) && (TS.m_iter%rom_interface.samplingFrequency() == 0)) {
        rom_interface.takeSample( sim, TS.m_waqt );
      }
#endif

      /* Call pre-step function */
      TimePreStep (&TS);
#ifdef compute_rhs_operators
      /* compute and write (to file) matrix operators representing the right-hand side */
//      if (((TS.m_iter+1)%solver->m_file_op_iter == 0) || (!TS.m_iter))
//        { ComputeRHSOperators(solver,mpi,TS.m_waqt);
#endif

      /* Step in time */
      TimeStep (&TS);

      /* Call post-step function */
      TimePostStep (&TS);

      ti_runtime += TS.m_iter_wctime;

      /* Print information to screen */
      TimePrintStep(&TS);

      /* Write intermediate solution to file */
      if (      ((TS.m_iter+1)%sim[0].solver.m_file_op_iter == 0)
            &&  ((TS.m_iter+1) < TS.m_n_iter) ) {
        for (int ns = 0; ns < a_nsims; ns++) {
          if (sim[ns].solver.PhysicsOutput) {
            sim[ns].solver.PhysicsOutput( &(sim[ns].solver),
                                          &(sim[ns].mpi),
                                          TS.m_waqt );
          }
        }
        OutputSolution(sim, a_nsims, TS.m_waqt);
#ifdef with_librom
        op_times_arr.push_back(TS.m_waqt);
#endif
      }

    }

    double t_final = TS.m_waqt;
    TimeCleanup(&TS);

    if (!a_rank) {
      printf( "Completed time integration (Final time: %f), total wctime: %f (seconds).\n",
              t_final, ti_runtime );
      if (a_nsims > 1) printf("\n");
    }

    /* calculate error if exact solution has been provided */
    for (int ns = 0; ns < a_nsims; ns++) {
      CalculateError(&(sim[ns].solver),
                     &(sim[ns].mpi) );
    }

    /* write a final solution file */
    for (int ns = 0; ns < a_nsims; ns++) {
      if (sim[ns].solver.PhysicsOutput) {
        sim[ns].solver.PhysicsOutput( &(sim[ns].solver),
                                      &(sim[ns].mpi),
                                      t_final );
      }
    }
    OutputSolution(sim, a_nsims, t_final);

#ifdef with_librom
    op_times_arr.push_back(TS.m_waqt);

    for (int ns = 0; ns < a_nsims; ns++) {
      ResetFilenameIndex( sim[ns].solver.m_filename_index,
                          sim[ns].solver.m_index_length );
    }

    if (rom_interface.mode() == _ROM_MODE_TRAIN_) {

      rom_interface.train();
      if (!a_rank) printf("libROM: total training wallclock time: %f (seconds).\n",
                        rom_interface.trainWallclockTime() );

      double total_rom_predict_time = 0;
      for (int iter = 0; iter < op_times_arr.size(); iter++) {

        double waqt = op_times_arr[iter];

        rom_interface.predict(sim, waqt);
        if (!a_rank) printf(  "libROM: Predicted solution at time %1.4e using ROM, wallclock time: %f.\n",
                            waqt, rom_interface.predictWallclockTime() );
        total_rom_predict_time += rom_interface.predictWallclockTime();

        /* calculate diff between ROM and PDE solutions */
        if (iter == (op_times_arr.size()-1)) {
          if (!a_rank) printf("libROM:   Calculating diff between PDE and ROM solutions.\n");
          for (int ns = 0; ns < a_nsims; ns++) {
            CalculateROMDiff(  &(sim[ns].solver),
                               &(sim[ns].mpi) );
          }
        }
        /* write the ROM solution to file */
        OutputROMSolution(sim, a_nsims,waqt);

      }

      if (!a_rank) {
        printf( "libROM: total prediction/query wallclock time: %f (seconds).\n",
                total_rom_predict_time );
      }

      rom_interface.saveROM();

    } else {

      for (int ns = 0; ns < a_nsims; ns++) {
        sim[ns].solver.m_rom_diff_norms[0]
          = sim[ns].solver.m_rom_diff_norms[1]
          = sim[ns].solver.m_rom_diff_norms[2]
          = -1;
      }

    }

  } else if (rom_mode == _ROM_MODE_PREDICT_) {

    for (int ns = 0; ns < a_nsims; ns++) {
      sim[ns].solver.m_rom_diff_norms[0]
        = sim[ns].solver.m_rom_diff_norms[1]
        = sim[ns].solver.m_rom_diff_norms[2]
        = -1;
      strcpy(sim[ns].solver.m_conservation_check,"no");
    }

    rom_interface.loadROM();
    rom_interface.projectInitialSolution(sim);

    {
      int start_iter = sim[0].solver.m_restart_iter;
      int n_iter = sim[0].solver.m_n_iter;
      double dt = sim[0].solver.m_dt;

      double cur_time = start_iter * dt;
      op_times_arr.push_back(cur_time);

      for (int iter = start_iter; iter < n_iter; iter++) {
        cur_time += dt;
        if (    ( (iter+1)%sim[0].solver.m_file_op_iter == 0)
            &&  ( (iter+1) < n_iter) ) {
          op_times_arr.push_back(cur_time);
        }
      }

      double t_final = n_iter*dt;
      op_times_arr.push_back(t_final);
    }

    double total_rom_predict_time = 0;
    for (int iter = 0; iter < op_times_arr.size(); iter++) {

      double waqt = op_times_arr[iter];

      rom_interface.predict(sim, waqt);
      if (!a_rank) printf(  "libROM: Predicted solution at time %1.4e using ROM, wallclock time: %f.\n",
                          waqt, rom_interface.predictWallclockTime() );
      total_rom_predict_time += rom_interface.predictWallclockTime();

      /* write the solution to file */
      for (int ns = 0; ns < a_nsims; ns++) {
        if (sim[ns].solver.PhysicsOutput) {
          sim[ns].solver.PhysicsOutput( &(sim[ns].solver),
                                        &(sim[ns].mpi),
                                        waqt );
        }
      }
      OutputSolution(sim, a_nsims, waqt);

    }

    /* calculate error if exact solution has been provided */
    for (int ns = 0; ns < a_nsims; ns++) {
      CalculateError(&(sim[ns].solver),
                     &(sim[ns].mpi) );
    }

    if (!a_rank) {
      printf( "libROM: total prediction/query wallclock time: %f (seconds).\n",
              total_rom_predict_time );
    }

  }
#endif

  return 0;
}
