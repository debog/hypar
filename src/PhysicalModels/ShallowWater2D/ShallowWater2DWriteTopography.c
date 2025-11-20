/*! @file ShallowWater2DWriteTopography.c
    @author Debojyoti Ghosh
    @brief Function to write out the topography
*/
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <common.h>
#include <arrayfunctions.h>
#include <io.h>
#include <plotfuncs.h>
#include <mpivars.h>
#include <hypar.h>
#include <physicalmodels/shallowwater2d.h>

/*! Write out the topography data to file */
int ShallowWater2DWriteTopography(  void*   s,  /*!< Solver object of type #HyPar */
                                    void*   m,  /*!< MPI object of type #MPIVariables */
                                    double  a_t /*!< Current simulation time */ )
{
  HyPar           *solver = (HyPar*)          s;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  ShallowWater2D  *params = (ShallowWater2D*) solver->m_physics;
  _DECLARE_IERR_;

  if (params->m_topo_flag) {

    char fname_root[_MAX_STRING_SIZE_] = "topography";
    if (solver->m_nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(solver->m_my_idx, index, (int)log10(solver->m_nsims)+1);
      strcat(fname_root, "_");
      strcat(fname_root, index);
      strcat(fname_root, "_");
    }

    WriteArray(  solver->m_ndims,
                 1,
                 solver->m_dim_global,
                 solver->m_dim_local,
                 solver->m_ghosts,
                 solver->m_x,
                 params->m_b,
                 solver,
                 mpi,
                 fname_root );

    if (!strcmp(solver->m_plot_solution, "yes")) {
      PlotArray( solver->m_ndims,
                 1,
                 solver->m_dim_global,
                 solver->m_dim_local,
                 solver->m_ghosts,
                 solver->m_x,
                 params->m_b,
                 a_t,
                 solver,
                 mpi,
                 fname_root );
    }

  }

  return(0);
}
