#ifdef with_librom

/*! @file OutputROMSolution.c
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

/* Function declarations */
void IncrementFilenameIndex(char*,int);

/*! Write out the ROM solution to file */
int OutputROMSolution(  void  *s,   /*!< Array of simulation objects of type #SimulationObject */
                        int   nsims /*!< Number of simulation objects */ )
{
  SimulationObject* simobj = (SimulationObject*) s;
  int ns;

  for (ns = 0; ns < nsims; ns++) {

    HyPar*        solver = &(simobj[ns].solver);
    MPIVariables* mpi    = &(simobj[ns].mpi);
    
    /* if WriteOutput() is NULL, then return */
    if (!solver->WriteOutput) continue;
  
    char fname_root[_MAX_STRING_SIZE_];
    strcpy(fname_root, solver->op_rom_fname_root);

    if (nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(ns, index, (int)log10(nsims)+1);
      strcat(fname_root, "_");
      strcat(fname_root, index);
    }
  
    WriteArray( solver->ndims,
                solver->nvars,
                solver->dim_global,
                solver->dim_local,
                solver->ghosts,
                solver->x,
                solver->u_rom_predicted,
                solver,
                mpi,
                fname_root );

    /* increment the index string, if required */
    if ((!strcmp(solver->output_mode,"serial")) && (!strcmp(solver->op_overwrite,"no"))) {
      IncrementFilenameIndex(solver->filename_index,solver->index_length);
    }

  }
  
  return 0;
}

#endif
