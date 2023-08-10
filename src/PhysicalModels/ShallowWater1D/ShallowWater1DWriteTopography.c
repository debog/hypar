/*! @file ShallowWater1DWriteTopography.c
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
#include <physicalmodels/shallowwater1d.h>

/*! Write out the topography data to file */
int ShallowWater1DWriteTopography(  void*   s,    /*!< Solver object of type #HyPar */
                                    void*   m,    /*!< MPI object of type #MPIVariables */
                                    double  a_t   /*!< Current simulation time */ )
{
  HyPar           *solver = (HyPar*)          s;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  ShallowWater1D  *params = (ShallowWater1D*) solver->physics;
  _DECLARE_IERR_;

  if (params->topo_flag) {

    char fname_root[_MAX_STRING_SIZE_] = "topography";
    if (solver->nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(solver->my_idx, index, (int)log10(solver->nsims)+1);
      strcat(fname_root, "_");
      strcat(fname_root, index);
      strcat(fname_root, "_");
    }

    WriteArray(  solver->ndims,
                 1,
                 solver->dim_global,
                 solver->dim_local,
                 solver->ghosts,
                 solver->x,
                 params->b,
                 solver,
                 mpi,
                 fname_root );

    if (!strcmp(solver->plot_solution, "yes")) {
      PlotArray( solver->ndims,
                 1,
                 solver->dim_global,
                 solver->dim_local,
                 solver->ghosts,
                 solver->x,
                 params->b,
                 a_t,
                 solver,mpi,
                 fname_root );
    }


  }

  return(0);
}
