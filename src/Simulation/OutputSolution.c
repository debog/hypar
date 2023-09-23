/*! @file OutputSolution.c
    @author Debojyoti Ghosh
    @brief Write out the solution to file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <common.h>
#include <arrayfunctions.h>
#include <io.h>
#include <timeintegration.h>
#include <mpivars.h>
#include <simulation_object.h>

#if defined(HAVE_CUDA)
# include <arrayfunctions_gpu.h>
#endif

/* Function declarations */
void IncrementFilenameIndex(char*,int);

/*! Write out the solution to file */
int OutputSolution( void  *s,   /*!< Array of simulation objects of type #SimulationObject */
                    int   nsims /*!< Number of simulation objects */
                  )
{
  SimulationObject* simobj = (SimulationObject*) s;
  int ns;
  _DECLARE_IERR_;

  for (ns = 0; ns < nsims; ns++) {

    HyPar*        solver = &(simobj[ns].solver);
    MPIVariables* mpi    = &(simobj[ns].mpi);

    /* if WriteOutput() is NULL, then return */
    if (!solver->WriteOutput) continue;

    /* time integration module may have auxiliary arrays to write out, so get them */
    int NSolutions = 0;
    IERR TimeGetAuxSolutions(&NSolutions,NULL,solver,-1,ns); CHECKERR(ierr);
    if (NSolutions > 10) NSolutions = 10;

    int  nu;
    char fname_root[_MAX_STRING_SIZE_];
    char aux_fname_root[_MAX_STRING_SIZE_];
    strcpy(fname_root, solver->op_fname_root);
    strcpy(aux_fname_root, solver->aux_op_fname_root);

    if (nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(ns, index, (int)log10(nsims)+1);
      strcat(fname_root, "_");
      strcat(fname_root, index);
      strcat(aux_fname_root, "_");
      strcat(aux_fname_root, index);
    }

    for (nu=0; nu<NSolutions; nu++) {

      double  *uaux = NULL;
      IERR TimeGetAuxSolutions(&NSolutions,&uaux,solver,nu,ns); CHECKERR(ierr);

      IERR WriteArray(  solver->ndims,
                        solver->nvars,
                        solver->dim_global,
                        solver->dim_local,
                        solver->ghosts,
                        solver->x,
                        uaux,
                        solver,
                        mpi,
                        aux_fname_root ); CHECKERR(ierr);

      aux_fname_root[2]++;
    }

#if defined(HAVE_CUDA)
    if (solver->use_gpu) {
      /* Copy values from GPU memory to CPU memory for writing. */
      gpuMemcpy(solver->x, solver->gpu_x, sizeof(double)*solver->size_x, gpuMemcpyDeviceToHost);

#ifdef CUDA_VAR_ORDERDING_AOS
      gpuMemcpy(  solver->u, 
                  solver->gpu_u, 
                  sizeof(double)*solver->ndof_cells_wghosts, 
                  gpuMemcpyDeviceToHost );
#else
      double *h_u = (double *) malloc(sizeof(double)*solver->ndof_cells_wghosts);
      gpuMemcpy(h_u, solver->gpu_u, sizeof(double)*solver->ndof_cells_wghosts, gpuMemcpyDeviceToHost);
      for (int i=0; i<solver->npoints_local_wghosts; i++) {
        for (int v=0; v<solver->nvars; v++) {
          solver->u[i*solver->nvars+v] = h_u[i+v*solver->npoints_local_wghosts];
        }
      }
      free(h_u);
#endif
    }
#endif

    IERR WriteArray(  solver->ndims,
                      solver->nvars,
                      solver->dim_global,
                      solver->dim_local,
                      solver->ghosts,
                      solver->x,
                      solver->u,
                      solver,
                      mpi,
                      fname_root ); CHECKERR(ierr);

    /* increment the index string, if required */
    if ((!strcmp(solver->output_mode,"serial")) && (!strcmp(solver->op_overwrite,"no"))) {
      IncrementFilenameIndex(solver->filename_index,solver->index_length);
    }

  }

  return(0);
}