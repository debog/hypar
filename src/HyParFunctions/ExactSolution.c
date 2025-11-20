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
  for (int i=0; i<solver->m_ndims; i++) {
    if (solver->m_dim_global[i] != solver->m_dim_global_ex[i]) same_size = 0;
  }

  if (same_size) {
    ReadArray( solver->m_ndims,
               solver->m_nvars,
               solver->m_dim_global,
               solver->m_dim_local,
               solver->m_ghosts,
               solver,
               m,
               NULL,
               uex,
               fname,
               flag);
  } else {
    ReadArraywInterp( solver->m_ndims,
                      solver->m_nvars,
                      solver->m_dim_global,
                      solver->m_dim_local,
                      solver->m_dim_global_ex,
                      solver->m_ghosts,
                      solver,
                      m,
                      NULL,
                      uex,
                      fname,
                      flag);
  }

  return(0);
}
