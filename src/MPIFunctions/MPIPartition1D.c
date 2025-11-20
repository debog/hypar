/*! @file MPIPartition1D.c
    @brief Compute the local size
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <mpivars.h>

/*!
  Given a 1D array of a given global size \a nglobal, and the total number
  of MPI ranks \a nproc on which it will be partitioned, this function
  computes the size of the local part of the 1D array on \a rank.
*/
int MPIPartition1D(
                    int a_nglobal,  /*!< Global size */
                    int a_nproc,    /*!< Total number of ranks */
                    int a_rank      /*!< Rank */
                  )
{
  int nlocal;
  if (a_nglobal%a_nproc == 0) nlocal = a_nglobal/a_nproc;
  else {
    if (a_rank == a_nproc-1)  nlocal = a_nglobal/a_nproc + a_nglobal%a_nproc;
    else                  nlocal = a_nglobal/a_nproc;
  }
  return(nlocal);
}
