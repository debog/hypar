/*! @file ApplyIBConditions.c
    @brief Apply immersed boundary conditions
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <mpivars.h>
#include <hypar.h>

/*!
  Call the physics-specific function that applies the immersed boundary conditions
  on the immersed boundary points.
*/
int ApplyIBConditions(void    *a_s, /*!< Object of type #HyPar containing solver-related variables */
                      void    *a_m, /*!< Object of type #MPIVariables containing MPI-related variables */
                      double  *a_x, /*!< The solution vector to apply immersed BCs on */
                      double  a_waqt /*!< Current simulation time */
                     )
{
  HyPar           *solver   = (HyPar*)          a_s;
  MPIVariables    *mpi      = (MPIVariables*)   a_m;

  /* Apply immersed boundary conditions, if applicable */
#if defined(HAVE_CUDA)
  if (solver->m_use_gpu) {
    if (solver->m_flag_ib) {
      fprintf(stderr, "ERROR: immersed boundaries have not yet been implemented on GPU.\n");
    }
  } else {
#endif
    if (solver->m_flag_ib) solver->IBFunction(solver,mpi,a_x,a_waqt);
#if defined(HAVE_CUDA)
  }
#endif

  return 0;
}
