/*! @file CompactSchemeInitialize.c
    @brief Initializes the compact schemes
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/*!
  This function initializes the compact finite-difference methods: allocates the arrays
  to store the tridiagonal system.
*/
int CompactSchemeInitialize(
                              void *s,      /*!< Solver object of type #HyPar */
                              void *m,      /*!< MPI object of type #MPIVariables */
                              char *type    /*!< Type of interpolation */
                           )
{
  HyPar         *solver = (HyPar*) s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  CompactScheme *compact   = (CompactScheme*) solver->m_compact;

  int nvars = solver->m_nvars;
  int ndims = solver->m_ndims;

  int size = 1, d;
  for (d=0; d<solver->m_ndims; d++) size *= (solver->m_dim_local[d]+1);
  size *= solver->m_nvars;
  if (!strcmp(solver->m_interp_type,_CHARACTERISTIC_)) size *= solver->m_nvars;

  compact->m_A = (double*) calloc (size, sizeof(double));
  compact->m_B = (double*) calloc (size, sizeof(double));
  compact->m_C = (double*) calloc (size, sizeof(double));
  compact->m_R = (double*) calloc (size, sizeof(double));

  compact->m_sendbuf = (double*) calloc (size, sizeof(double));
  compact->m_recvbuf = (double*) calloc (size, sizeof(double));

  return(0);
}
