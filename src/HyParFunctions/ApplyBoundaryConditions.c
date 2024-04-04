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
 * (#HyPar::boundary[#HyPar::nBoundaryZones]) and calls the corresponding boundary
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
  DomainBoundary  *boundary = (DomainBoundary*) solver->boundary;
  MPIVariables    *mpi      = (MPIVariables*)   m;
  int             nb        = solver->nBoundaryZones;

  int* dim_local;
#if defined(HAVE_CUDA)
  if (solver->use_gpu) {
    dim_local = solver->gpu_dim_local;
  } else {
#endif
    dim_local = solver->dim_local;
#if defined(HAVE_CUDA)
  }
#endif

  /* Apply domain boundary conditions to x */
  int n;
  for (n = 0; n < nb; n++) {
    boundary[n].BCFunctionU(&boundary[n],mpi,solver->ndims,solver->nvars,
                            dim_local,solver->ghosts,x,waqt);
  }

  return(0);
}
