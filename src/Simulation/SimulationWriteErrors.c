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
#include <simulation.h>

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
    for (n = 0; n < nsims; n++) {

      char  err_fname[_MAX_STRING_SIZE_],
            cons_fname[_MAX_STRING_SIZE_],
            fc_fname[_MAX_STRING_SIZE_];

      strcpy(err_fname,"errors");
      strcpy(cons_fname,"conservation");
      strcpy(fc_fname,"function_counts");

      if (nsims > 1) {

        strcat(err_fname,"_");
        strcat(cons_fname,"_");
        strcat(fc_fname,"_");

        char index[_MAX_STRING_SIZE_];
        GetStringFromInteger(n, index,4);

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
      for (d=0; d<sim[n].solver.ndims; d++) fprintf(out,"%4d ",sim[n].solver.dim_global[d]);
      for (d=0; d<sim[n].solver.ndims; d++) fprintf(out,"%4d ",sim[n].mpi.iproc[d]);
      fprintf(out,"%1.16E  ",sim[n].solver.dt);
      fprintf(out,"%1.16E %1.16E %1.16E   ",sim[n].solver.error[0],sim[n].solver.error[1],sim[n].solver.error[2]);
      fprintf(out,"%1.16E %1.16E\n",solver_runtime,main_runtime);
      fclose(out);
      /* write out conservation errors to file */
      out = fopen(cons_fname,"w");
      for (d=0; d<sim[n].solver.ndims; d++) fprintf(out,"%4d ",sim[n].solver.dim_global[d]);
      for (d=0; d<sim[n].solver.ndims; d++) fprintf(out,"%4d ",sim[n].mpi.iproc[d]);
      fprintf(out,"%1.16E  ",sim[n].solver.dt);
      for (d=0; d<sim[n].solver.nvars; d++) fprintf(out,"%1.16E ",sim[n].solver.ConservationError[d]);
      fprintf(out,"\n");
      fclose(out);
      /* write out function call counts to file */
      out = fopen(fc_fname,"w");
      fprintf(out,"%d\n",sim[n].solver.n_iter);
      fprintf(out,"%d\n",sim[n].solver.count_hyp);
      fprintf(out,"%d\n",sim[n].solver.count_par);
      fprintf(out,"%d\n",sim[n].solver.count_sou);
#ifdef with_petsc
      fprintf(out,"%d\n",sim[n].solver.count_RHSFunction);
      fprintf(out,"%d\n",sim[n].solver.count_IFunction);
      fprintf(out,"%d\n",sim[n].solver.count_IJacobian);
      fprintf(out,"%d\n",sim[n].solver.count_IJacFunction);
#endif
      fclose(out);

      /* print solution errors, conservation errors, and wall times to screen */
      printf("Computed errors for simuation %3d:\n", n);
      printf("  L1         Error           : %1.16E\n",sim[n].solver.error[0]);
      printf("  L2         Error           : %1.16E\n",sim[n].solver.error[1]);
      printf("  Linfinity  Error           : %1.16E\n",sim[n].solver.error[2]);
      printf("Conservation Errors:\n");
      for (d=0; d<sim[n].solver.nvars; d++) printf("\t%1.16E\n",sim[n].solver.ConservationError[d]);
      printf("Solver runtime (in seconds): %1.16E\n",solver_runtime);
      printf("Total  runtime (in seconds): %1.16E\n",main_runtime);

    }
  }

  return;
}
