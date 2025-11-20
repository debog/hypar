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
int OutputROMSolution(  void*   a_s,      /*!< Array of simulation objects of type #SimulationObject */
                        int     a_nsims,  /*!< Number of simulation objects */
                        double  a_time  /*!< Current simulation time */)
{
  SimulationObject* simobj = (SimulationObject*) a_s;
  int ns;

  for (ns = 0; ns < a_nsims; ns++) {

    HyPar*        solver = &(simobj[ns].solver);
    MPIVariables* mpi    = &(simobj[ns].mpi);

    if ((!solver->WriteOutput) && (strcmp(solver->m_plot_solution,"yes"))) continue;

    char fname_root[_MAX_STRING_SIZE_];
    strcpy(fname_root, solver->m_op_rom_fname_root);

    if (a_nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(ns, index, (int)log10(a_nsims)+1);
      strcat(fname_root, "_");
      strcat(fname_root, index);
    }

    WriteArray( solver->m_ndims,
                solver->m_nvars,
                solver->m_dim_global,
                solver->m_dim_local,
                solver->m_ghosts,
                solver->m_x,
                solver->m_u_rom_predicted,
                solver,
                mpi,
                fname_root );

    if (!strcmp(solver->m_plot_solution, "yes")) {
      PlotArray(   solver->m_ndims,
                   solver->m_nvars,
                   solver->m_dim_global,
                   solver->m_dim_local,
                   solver->m_ghosts,
                   solver->m_x,
                   solver->m_u,
                   a_time,
                   solver,
                   mpi,
                   fname_root );
    }

    /* increment the index string, if required */
    if ((!strcmp(solver->m_output_mode,"serial")) && (!strcmp(solver->m_op_overwrite,"no"))) {
      IncrementFilenameIndex(solver->m_filename_index,solver->m_index_length);
    }

  }

  return 0;
}

#endif
