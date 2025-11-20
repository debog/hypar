/*! @file MPIPartitionArray1D.c
    @brief Partition an essentially 1D array
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

/*!
  Partitions the contents of a global 1D array on the root rank (rank 0) to local arrays
  on all the MPI ranks. See documentation of MPIExchangeBoundaries1D() on what a "1D
  array" is in the context of a multidimensional simulation. The 1D array must be the same
  along spatial dimensions normal to the one it represents.

  Notes:
  + The global array must not have ghost points.
  + The global array must be preallocated on only rank 0. On other ranks, it must be NULL.
  + Since this function deals with a 1D array, more than one rank may be receiving the same
    piece of data from rank 0 (i.e. if there are more than one MPI rank along the dimensions
    normal to one corresponding to \a x ).
*/
int MPIPartitionArray1D(
                          void    *a_m,       /*!< MPI object of type #MPIVariables */
                          double  *a_xg,      /*!< Global 1D array (must be preallocated) without ghost points */
                          double  *a_x,       /*!< Local 1D array to be gathered */
                          int     a_istart,   /*!< Starting index (global) of this rank's portion of the array */
                          int     a_iend,     /*!< Ending index (global) of this rank's portion of the array + 1 */
                          int     a_N_local,  /*!< Local size of the array */
                          int     a_ghosts    /*!< Number of ghost points */
                        )
{
  MPIVariables *mpi = (MPIVariables*) a_m;
  int          ierr = 0;
#ifndef serial
  MPI_Status   status;
#endif

  /* a_xg should be non-null only on root */
  if (mpi->m_rank && a_xg) {
    fprintf(stderr,"Error in MPIPartitionArray1D(): global array exists on non-root processors (rank %d).\n",
            mpi->m_rank);
    ierr = 1;
  }
  if ((!mpi->m_rank) && (!a_xg)) {
    fprintf(stderr,"Error in MPIPartitionArray1D(): global array is not allocated on root processor.\n");
    ierr = 1;
  }

  if (!mpi->m_rank) {
    int proc;
    for (proc = 0; proc < mpi->m_nproc; proc++) {
      /* Find out the domain limits for each process */
      int is,ie;
      if (proc) {
#ifndef serial
        MPI_Recv(&is,1,MPI_INT,proc,1442,mpi->m_world,&status);
        MPI_Recv(&ie,1,MPI_INT,proc,1443,mpi->m_world,&status);
#endif
      } else {
        is = a_istart;
        ie = a_iend;
      }
      int size = ie - is;
      double *buffer = (double*) calloc (size, sizeof(double));
      _ArrayCopy1D_((a_xg+is),buffer,size);
      if (proc) {
#ifndef serial
        MPI_Send(buffer,size,MPI_DOUBLE,proc,1539,mpi->m_world);
#endif
      } else _ArrayCopy1D_(buffer,a_x,a_N_local);
      free(buffer);
    }

  } else {
#ifndef serial
    /* Meanwhile, on other processes */
    /* send local start and end indices to root */
    MPI_Send(&a_istart,1,MPI_INT,0,1442,mpi->m_world);
    MPI_Send(&a_iend  ,1,MPI_INT,0,1443,mpi->m_world);
    double *buffer = (double*) calloc (a_N_local, sizeof(buffer));
    MPI_Recv(buffer,a_N_local,MPI_DOUBLE,0,1539,mpi->m_world,&status);
    _ArrayCopy1D_(buffer,a_x,a_N_local);
    free(buffer);
#endif
  }
  return(ierr);
}
