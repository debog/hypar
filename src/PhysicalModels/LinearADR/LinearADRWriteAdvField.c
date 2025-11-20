/*! @file LinearADRWriteAdvField.c
    @author Debojyoti Ghosh
    @brief Write out the advection field to file
*/
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <common.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <io.h>
#include <plotfuncs.h>
#include <mpivars.h>
#include <hypar.h>
#include <physicalmodels/linearadr.h>

/*! Write out the advection field to file */
int LinearADRWriteAdvField( void*   a_s,    /*!< Solver object of type #HyPar */
                            void*   a_m,    /*!< MPI object of type #MPIVariables */
                            double  a_t   /*!< Current simulation time */ )
{
  HyPar         *solver = (HyPar*)          a_s;
  MPIVariables  *mpi    = (MPIVariables*)   a_m;
  LinearADR     *param  = (LinearADR*) solver->m_physics;

  if (param->m_constant_advection == 0) {

    char fname_root[_MAX_STRING_SIZE_] = "advection_field";
    if (solver->m_nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(solver->m_my_idx, index, (int)log10(solver->m_nsims)+1);
      strcat(fname_root, "_");
      strcat(fname_root, index);
    }

    int adv_nvar = solver->m_ndims * solver->m_nvars;
    WriteArray(  solver->m_ndims,
                 adv_nvar,
                 solver->m_dim_global,
                 solver->m_dim_local,
                 solver->m_ghosts,
                 solver->m_x,
                 param->m_a,
                 solver,mpi,
                 fname_root );

    if (!strcmp(solver->m_plot_solution, "yes")) {
      PlotArray( solver->m_ndims,
                 adv_nvar,
                 solver->m_dim_global,
                 solver->m_dim_local,
                 solver->m_ghosts,
                 solver->m_x,
                 param->m_a,
                 a_t,
                 solver,mpi,
                 fname_root );
    }

  }

  return(0);
}
