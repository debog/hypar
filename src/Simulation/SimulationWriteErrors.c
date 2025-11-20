/*! @file SimulationWriteErrors.c
    @brief Write errors for each simulation
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <basic.h>
#include <common.h>
#include <mpivars.h>
#include <simulation_object.h>

/*! Writes out the errors and other data for each simulation.
*/

void SimWriteErrors(void  *s,               /*!< Array of simulations of type #SimulationObject */
                    int   nsims,            /*!< Number of simulations */
                    int   rank,             /*!< MPI rank of this process */
                    double solver_runtime,  /*!< Measured runtime of solver */
                    double main_runtime     /*!< Measured total runtime */
                   )
{
  SimulationObject* sim = (SimulationObject*) s;
  int n;

  if (!rank) {

    if (nsims > 1) printf("\n");

    for (n = 0; n < nsims; n++) {

      char err_fname[_MAX_STRING_SIZE_],
           cons_fname[_MAX_STRING_SIZE_],
           fc_fname[_MAX_STRING_SIZE_];
      strcpy(err_fname,"errors");
      strcpy(cons_fname,"conservation");
      strcpy(fc_fname,"function_counts");
#ifdef with_librom
      char rom_diff_fname[_MAX_STRING_SIZE_];
      strcpy(rom_diff_fname,"pde_rom_diff");
#endif


      if (nsims > 1) {

        strcat(err_fname,"_");
        strcat(cons_fname,"_");
        strcat(fc_fname,"_");
#ifdef with_librom
        strcat(rom_diff_fname,"_");
#endif

        char index[_MAX_STRING_SIZE_];
        GetStringFromInteger(n, index, (int)log10(nsims)+1);

        strcat(err_fname,index);
        strcat(cons_fname,index);
        strcat(fc_fname,index);
#ifdef with_librom
        strcat(rom_diff_fname,index);
#endif
      }

      strcat(err_fname,".dat");
      strcat(cons_fname,".dat");
      strcat(fc_fname,".dat");
#ifdef with_librom
      strcat(rom_diff_fname,".dat");
#endif

      FILE *out;
      /* write out solution errors and wall times to file */
      out = fopen(err_fname,"w");
      for (int d=0; d<sim[n].solver.m_ndims; d++) fprintf(out,"%4d ",sim[n].solver.m_dim_global[d]);
      for (int d=0; d<sim[n].solver.m_ndims; d++) fprintf(out,"%4d ",sim[n].mpi.m_iproc[d]);
      fprintf(out,"%1.16E  ",sim[n].solver.m_dt);
      fprintf(out,"%1.16E %1.16E %1.16E   ",sim[n].solver.m_error[0],sim[n].solver.m_error[1],sim[n].solver.m_error[2]);
      fprintf(out,"%1.16E %1.16E\n",solver_runtime,main_runtime);
      fclose(out);
      /* write out conservation errors to file */
      out = fopen(cons_fname,"w");
      for (int d=0; d<sim[n].solver.m_ndims; d++) fprintf(out,"%4d ",sim[n].solver.m_dim_global[d]);
      for (int d=0; d<sim[n].solver.m_ndims; d++) fprintf(out,"%4d ",sim[n].mpi.m_iproc[d]);
      fprintf(out,"%1.16E  ",sim[n].solver.m_dt);
      for (int d=0; d<sim[n].solver.m_nvars; d++) fprintf(out,"%1.16E ",sim[n].solver.m_conservation_error[d]);
      fprintf(out,"\n");
      fclose(out);
      /* write out function call counts to file */
      out = fopen(fc_fname,"w");
      fprintf(out,"%d\n",sim[n].solver.m_n_iter);
      fprintf(out,"%d\n",sim[n].solver.m_count_hyp);
      fprintf(out,"%d\n",sim[n].solver.m_count_par);
      fprintf(out,"%d\n",sim[n].solver.m_count_sou);
#ifdef with_petsc
      fprintf(out,"%d\n",sim[n].solver.m_count_rhs_function);
      fprintf(out,"%d\n",sim[n].solver.m_count_i_function);
      fprintf(out,"%d\n",sim[n].solver.m_count_i_jacobian);
      fprintf(out,"%d\n",sim[n].solver.m_count_i_jac_function);
#endif
      fclose(out);
#ifdef with_librom
      /* write out solution errors and wall times to file */
      if (sim[n].solver.m_rom_diff_norms[0] >= 0) {
        out = fopen(rom_diff_fname,"w");
        for (int d=0; d<sim[n].solver.m_ndims; d++) fprintf(out,"%4d ",sim[n].solver.m_dim_global[d]);
        for (int d=0; d<sim[n].solver.m_ndims; d++) fprintf(out,"%4d ",sim[n].mpi.m_iproc[d]);
        fprintf(out,"%1.16E  ",sim[n].solver.m_dt);
        fprintf(out,"%1.16E %1.16E %1.16E   ",sim[n].solver.m_rom_diff_norms[0],sim[n].solver.m_rom_diff_norms[1],sim[n].solver.m_rom_diff_norms[2]);
        fprintf(out,"%1.16E %1.16E\n",solver_runtime,main_runtime);
        fclose(out);
      }
#endif

      /* print solution errors, conservation errors, and wall times to screen */
      if (sim[n].solver.m_error[0] >= 0) {
        printf("Computed errors for domain %d:\n", n);
        printf("  L1         Error           : %1.16E\n",sim[n].solver.m_error[0]);
        printf("  L2         Error           : %1.16E\n",sim[n].solver.m_error[1]);
        printf("  Linfinity  Error           : %1.16E\n",sim[n].solver.m_error[2]);
      }
      if (!strcmp(sim[n].solver.m_conservation_check,"yes")) {
        printf("Conservation Errors:\n");
        for (int d=0; d<sim[n].solver.m_nvars; d++) printf("\t%1.16E\n",sim[n].solver.m_conservation_error[d]);
      }
#ifdef with_librom
      if (sim[n].solver.m_rom_diff_norms[0] >= 0) {
        printf("Norms of the diff between ROM and PDE solutions for domain %d:\n", n);
        printf("  L1         Norm            : %1.16E\n",sim[n].solver.m_rom_diff_norms[0]);
        printf("  L2         Norm            : %1.16E\n",sim[n].solver.m_rom_diff_norms[1]);
        printf("  Linfinity  Norm            : %1.16E\n",sim[n].solver.m_rom_diff_norms[2]);
      }
#endif

    }

    printf("Solver runtime (in seconds): %1.16E\n",solver_runtime);
    printf("Total  runtime (in seconds): %1.16E\n",main_runtime);
    if (nsims > 1) printf("\n");

  }

  return;
}
