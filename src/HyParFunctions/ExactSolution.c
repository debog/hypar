/*! @file ExactSolution.c
    @author Debojyoti Ghosh
    @brief Read in exact solution, if available.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <io.h>
#include <mpivars.h>
#include <hypar.h>

/*! Read in the exact solution, if available. */
int ExactSolution(
                    void   *s,     /*!< Solver object of type #HyPar */
                    void   *m,     /*!< MPI object of type #MPIVariables */
                    double *uex,   /*!< Array to hold the exact solution, if available */
                    char   *fname, /*!< Filename root from which to read exact solution */
                    int    *flag   /*!< Flag to indicate if exact solution was available */
                 )
{
  HyPar  *solver = (HyPar*) s;

  int same_size = 1;
  for (int i=0; i<solver->ndims; i++) {
    if (solver->dim_global[i] != solver->dim_global_ex[i]) same_size = 0;
  }

  if (same_size) {
    ReadArray( solver->ndims,
               solver->nvars,
               solver->dim_global,
               solver->dim_local,
               solver->ghosts,
               solver,
               m,
               NULL,
               uex,
               fname,
               flag);
  } else {
    ReadArraywInterp( solver->ndims,
                      solver->nvars,
                      solver->dim_global,
                      solver->dim_local,
                      solver->dim_global_ex,
                      solver->ghosts,
                      solver,
                      m,
                      NULL,
                      uex,
                      fname,
                      flag);
  }

  return(0);
}
