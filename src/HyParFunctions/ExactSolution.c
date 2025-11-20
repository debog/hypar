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
                    void   *a_s,     /*!< Solver object of type #HyPar */
                    void   *a_m,     /*!< MPI object of type #MPIVariables */
                    double *a_uex,   /*!< Array to hold the exact solution, if available */
                    char   *a_fname, /*!< Filename root from which to read exact solution */
                    int    *a_flag   /*!< Flag to indicate if exact solution was available */
                 )
{
  HyPar  *solver = (HyPar*) a_s;

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
               a_m,
               NULL,
               a_uex,
               a_fname,
               a_flag);
  } else {
    ReadArraywInterp( solver->m_ndims,
                      solver->m_nvars,
                      solver->m_dim_global,
                      solver->m_dim_local,
                      solver->m_dim_global_ex,
                      solver->m_ghosts,
                      solver,
                      a_m,
                      NULL,
                      a_uex,
                      a_fname,
                      a_flag);
  }

  return(0);
}
