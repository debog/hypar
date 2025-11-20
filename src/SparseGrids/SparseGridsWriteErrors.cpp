/*! @file SparseGridsWriteErrors.cpp
    @brief Write errors for the simulation
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <string>
#include <basic.h>
#include <common_cpp.h>
#include <sparse_grids_simulation.h>


/*! Writes out the errors and other data for the sparse grids simulation. */
void SparseGridsSimulation::WriteErrors(double solver_runtime,  /*!< Measured runtime of solver */
                                        double main_runtime     /*!< Measured total runtime */
                                       )
{
  if (!m_rank) {

    /* Write sparse grids stuff, if asked for */
    if (m_print_sg_errors == 1) {
      for (int n = 0; n < m_nsims_sg; n++) {

        char  err_fname[_MAX_STRING_SIZE_],
              cons_fname[_MAX_STRING_SIZE_],
              fc_fname[_MAX_STRING_SIZE_];

        strcpy(err_fname, "errors");
        strcpy(cons_fname,"conservation");
        strcpy(fc_fname,  "function_counts");

        if (m_nsims_sg > 1) {

          strcat(err_fname,"_");
          strcat(cons_fname,"_");
          strcat(fc_fname,"_");

          char index[_MAX_STRING_SIZE_];
          GetStringFromInteger(n, index, (int)log10(m_nsims_sg)+1);

          strcat(err_fname,index);
          strcat(cons_fname,index);
          strcat(fc_fname,index);
        }

        strcat(err_fname,".dat");
        strcat(cons_fname,".dat");
        strcat(fc_fname,".dat");

        FILE *out;
        /* write out solution errors and wall times to file */
        int d;
        out = fopen(err_fname,"w");
        for (d=0; d<m_sims_sg[n].solver.m_ndims; d++) fprintf(out,"%4d ",m_sims_sg[n].solver.m_dim_global[d]);
        for (d=0; d<m_sims_sg[n].solver.m_ndims; d++) fprintf(out,"%4d ",m_sims_sg[n].mpi.m_iproc[d]);
        fprintf(out,"%1.16E  ",m_sims_sg[n].solver.m_dt);
        fprintf(out,"%1.16E %1.16E %1.16E   ",m_sims_sg[n].solver.m_error[0],m_sims_sg[n].solver.m_error[1],m_sims_sg[n].solver.m_error[2]);
        fprintf(out,"%1.16E %1.16E\n",solver_runtime,main_runtime);
        fclose(out);
        /* write out conservation errors to file */
        out = fopen(cons_fname,"w");
        for (d=0; d<m_sims_sg[n].solver.m_ndims; d++) fprintf(out,"%4d ",m_sims_sg[n].solver.m_dim_global[d]);
        for (d=0; d<m_sims_sg[n].solver.m_ndims; d++) fprintf(out,"%4d ",m_sims_sg[n].mpi.m_iproc[d]);
        fprintf(out,"%1.16E  ",m_sims_sg[n].solver.m_dt);
        for (d=0; d<m_sims_sg[n].solver.m_nvars; d++) fprintf(out,"%1.16E ",m_sims_sg[n].solver.m_conservation_error[d]);
        fprintf(out,"\n");
        fclose(out);
        /* write out function call counts to file */
        out = fopen(fc_fname,"w");
        fprintf(out,"%d\n",m_sims_sg[n].solver.m_n_iter);
        fprintf(out,"%d\n",m_sims_sg[n].solver.m_count_hyp);
        fprintf(out,"%d\n",m_sims_sg[n].solver.m_count_par);
        fprintf(out,"%d\n",m_sims_sg[n].solver.m_count_sou);
  #ifdef with_petsc
        fprintf(out,"%d\n",m_sims_sg[n].solver.m_count_rhs_function);
        fprintf(out,"%d\n",m_sims_sg[n].solver.m_count_i_function);
        fprintf(out,"%d\n",m_sims_sg[n].solver.m_count_i_jacobian);
        fprintf(out,"%d\n",m_sims_sg[n].solver.m_count_i_jac_function);
  #endif
        fclose(out);

        /* print solution errors, conservation errors, and wall times to screen */
        if (m_sims_sg[n].solver.m_error[0] >= 0) {
          printf("Computed errors for sparse grids domain %d:\n", n);
          printf("  L1   Error: %1.16E\n",m_sims_sg[n].solver.m_error[0]);
          printf("  L2   Error: %1.16E\n",m_sims_sg[n].solver.m_error[1]);
          printf("  Linf Error: %1.16E\n",m_sims_sg[n].solver.m_error[2]);
        }
        if (!strcmp(m_sims_sg[n].solver.m_conservation_check,"yes")) {
          printf("Conservation Errors:\n");
          for (d=0; d<m_sims_sg[n].solver.m_nvars; d++) printf("\t%1.16E\n",m_sims_sg[n].solver.m_conservation_error[d]);
          printf("\n");
        }

      }
    }

    /* First write stuff for the full grid solution */
    {
      char  err_fname[_MAX_STRING_SIZE_];
      strcpy(err_fname,"errors_fg");
      strcat(err_fname,".dat");

      FILE *out;
      /* write out solution errors and wall times to file */
      int d;
      out = fopen(err_fname,"w");
      for (d=0; d<m_sim_fg->solver.m_ndims; d++) fprintf(out,"%4d ",m_sim_fg->solver.m_dim_global[d]);
      for (d=0; d<m_sim_fg->solver.m_ndims; d++) fprintf(out,"%4d ",m_sim_fg->mpi.m_iproc[d]);
      fprintf(out,"%1.16E  ",m_sim_fg->solver.m_dt);
      fprintf(out,"%1.16E %1.16E %1.16E   ",m_sim_fg->solver.m_error[0],m_sim_fg->solver.m_error[1],m_sim_fg->solver.m_error[2]);
      fprintf(out,"%1.16E %1.16E\n",solver_runtime,main_runtime);
      fclose(out);

      /* print solution errors, conservation errors, and wall times to screen */
      if (m_sim_fg->solver.m_error[0] >= 0) {
        printf("Computed errors for full grid solution:\n");
        printf("  L1   Error: %1.16E\n",m_sim_fg->solver.m_error[0]);
        printf("  L2   Error: %1.16E\n",m_sim_fg->solver.m_error[1]);
        printf("  Linf Error: %1.16E\n",m_sim_fg->solver.m_error[2]);
      }
    }

    /* report wall times */
    printf("Solver runtime (in seconds): %1.16E\n",solver_runtime);
    printf("Total  runtime (in seconds): %1.16E\n",main_runtime);

  }

  return;
}
