#ifdef with_librom

/*! @file OutputROMSolution.cpp
    @author Debojyoti Ghosh
    @brief Write out the solution to file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <common_cpp.h>
#include <arrayfunctions.h>
#include <io_cpp.h>
#include <plotfuncs_cpp.h>
#include <mpivars_cpp.h>
#include <simulation_object.h>

/* Function declarations */
extern "C" void IncrementFilenameIndex(char*,int);

/*! Write out the ROM solution to file */
int OutputROMSolution(  void*   s,      /*!< Array of simulation objects of type #SimulationObject */
                        int     nsims,  /*!< Number of simulation objects */
                        double  a_time  /*!< Current simulation time */)
{
  SimulationObject* simobj = (SimulationObject*) s;
  int ns;

  for (ns = 0; ns < nsims; ns++) {

    HyPar*        solver = &(simobj[ns].solver);
    MPIVariables* mpi    = &(simobj[ns].mpi);

    if ((!solver->WriteOutput) && (strcmp(solver->plot_solution,"yes"))) continue;

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

    if (!strcmp(solver->plot_solution, "yes")) {
      PlotArray(   solver->ndims,
                   solver->nvars,
                   solver->dim_global,
                   solver->dim_local,
                   solver->ghosts,
                   solver->x,
                   solver->u,
                   a_time,
                   solver,
                   mpi,
                   fname_root );
    }

    /* increment the index string, if required */
    if ((!strcmp(solver->output_mode,"serial")) && (!strcmp(solver->op_overwrite,"no"))) {
      IncrementFilenameIndex(solver->filename_index,solver->index_length);
    }

  }

  return 0;
}

#endif
