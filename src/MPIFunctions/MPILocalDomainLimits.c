/*! @file MPILocalDomainLimits.c
    @brief Compute local domain limites of a MPI rank
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <basic.h>
#include <mpivars.h>

/*!
  Computes the local domain limites: Given a MPI rank \a p, this function will compute the starting and ending indices
  of the local domain on this rank. These indices are global.
  + The starting and ending indices are stored in preallocated integer arrays, whose elements are these indices in each
    spatial dimension.
  + The ending index is 1 + actual last index; thus the one can write the loop as (i=start; i<end; i++) and \b not
    (i=start; i<=end; i++).
*/
int MPILocalDomainLimits(
                          int   a_ndims,        /*!< Number of spatial dimensions */
                          int   a_p,            /*!< MPI rank */
                          void  *a_m,           /*!< MPI object of type #MPIVariables */
                          int   *a_dim_global,  /*!< Integer array with elements as global size in each spatial dimension */
                          int   *a_is,          /*!< Integer array whose elements will contain the starting index of the local domain on rank \a p */
                          int   *a_ie           /*!< Integer array whose elements will contain the ending index of the local domain on rank \a p */
                        )
{
  MPIVariables *mpi = (MPIVariables*) a_m;
  int          i;
  _DECLARE_IERR_;

  int ip[a_ndims];
  IERR MPIRanknD(a_ndims,a_p,mpi->m_iproc,ip); CHECKERR(ierr);

  for (i=0; i<a_ndims; i++) {
    int imax_local, isize, root = 0;
    imax_local = MPIPartition1D(a_dim_global[i],mpi->m_iproc[i],root );
    isize      = MPIPartition1D(a_dim_global[i],mpi->m_iproc[i],ip[i]);
    if (a_is)  a_is[i] = ip[i]*imax_local;
    if (a_ie)  a_ie[i] = ip[i]*imax_local + isize;
  }
  return(0);
}
