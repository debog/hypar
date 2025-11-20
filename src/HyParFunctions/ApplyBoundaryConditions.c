/*! @file ApplyBoundaryConditions.c
 *  @author Debojyoti Ghosh
 *  @brief Apply physical boundary conditions to domain.
 *
 *  Contains the function that applies the physical boundary conditions
 *  to each boundary zone.
*/

#include <stdio.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <boundaryconditions.h>
#include <hypar.h>

/*!
 * \brief Applies the boundary conditions specified for each boundary zone.
 *
 * The solver object (of type #HyPar) contains an oject of type #DomainBoundary
 * that contains all the boundary information (dimension, extent, face, type, etc).
 * This function iterates through each of the boundary zones
 * (#HyPar::m_boundary[#HyPar::m_n_boundary_zones]) and calls the corresponding boundary
 * condition function.
 * \n\n
 * The variable \a flag indicates if the array \a x is the solution, or a delta-solution
 * (from implicit time-integration methods).
*/
int ApplyBoundaryConditions(void    *s,     /*!< Object of type #HyPar containing solver-related variables */
                            void    *m,     /*!< Object of type #MPIVariables containing MPI-related variables */
                            double  *x,     /*!< The solution vector on which the boundary conditions are to be applied */
                            double  *xref,  /*!< Reference solution vector, if needed */
                            double  waqt    /*!< Current simulation time */
                           )
{
  HyPar           *solver   = (HyPar*)          s;
  DomainBoundary  *boundary = (DomainBoundary*) solver->m_boundary;
  MPIVariables    *mpi      = (MPIVariables*)   m;
  int             nb        = solver->m_n_boundary_zones;

  int* dim_local;
#if defined(HAVE_CUDA)
  if (solver->m_use_gpu) {
    dim_local = solver->m_gpu_dim_local;
  } else {
#endif
    dim_local = solver->m_dim_local;
#if defined(HAVE_CUDA)
  }
#endif

  /* Apply domain boundary conditions to x */
  int n;
  for (n = 0; n < nb; n++) {
    boundary[n].BCFunctionU(&boundary[n],mpi,solver->m_ndims,solver->m_nvars,
                            dim_local,solver->m_ghosts,x,waqt);
  }

  return(0);
}
