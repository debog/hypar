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
#include <mpivars.h>
#include <hypar.h>
#include <physicalmodels/linearadr.h>

/*! Write out the advection field to file */
int LinearADRWriteAdvField( void* s, /*!< Solver object of type #HyPar */
                            void* m  /*!< MPI object of type #MPIVariables */
                          )
{
  HyPar         *solver = (HyPar*)          s;
  MPIVariables  *mpi    = (MPIVariables*)   m;
  LinearADR     *param  = (LinearADR*) solver->physics;

  if (param->constant_advection == 0) {

    char fname_root[_MAX_STRING_SIZE_] = "advection_field";
    if (solver->nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(solver->my_idx, index, (int)log10(solver->nsims)+1);
      strcat(fname_root, "_");
      strcat(fname_root, index);
      strcat(fname_root, "_");
    }

    int adv_nvar = solver->ndims * solver->nvars;
    IERR WriteArray(  solver->ndims,
                      adv_nvar,
                      solver->dim_global,
                      solver->dim_local,
                      solver->ghosts,
                      solver->x,
                      param->a,
                      solver,mpi,
                      fname_root ); CHECKERR(ierr);

  }

  return(0);
}
