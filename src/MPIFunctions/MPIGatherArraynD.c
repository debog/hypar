/*! @file MPIGatherArraynD.c
    @brief Gather an n-dimensional array in to a global array
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

/*!
  Gathers the contents of an n-dimensional array, partitioned amongst the MPI ranks, in to a global
  array on rank 0. See documentation of MPIExchangeBoundariesnD() for how the n-dimensional array is
  stored in the memory as a single-index array.

  Notes:
  + The global array must have no ghost points.
  + The global array must be allocated only on rank 0. On other ranks, it must be NULL.
*/
int MPIGatherArraynD(
                      int     a_ndims,        /*!< Number of spatial dimensions */
                      void    *a_m,           /*!< MPI object of type #MPIVariables */
                      double  *a_xg,          /*!< Global array (preallocated) without ghost points */
                      double  *a_x,           /*!< Local array */
                      int     *a_dim_global,  /*!< Integer array with elements as global size along each spatial dimension */
                      int     *a_dim_local,   /*!< Integer array with elements as local size along each spatial dimension */
                      int     a_ghosts,       /*!< Number of ghost points */
                      int     a_nvars         /*!< Number of variables (vector components) */
                    )
{
  MPIVariables *mpi = (MPIVariables*) a_m;
  int          d, size;
  _DECLARE_IERR_;

  int is[a_ndims], ie[a_ndims], index[a_ndims], bounds[a_ndims];

  /* a_xg should be non-null only on root */
  if (mpi->m_rank && a_xg) {
    fprintf(stderr,"Error in MPIGatherArraynD(): global array exists on non-root processors (rank %d).\n",
            mpi->m_rank);
    return(1);
  }
  if ((!mpi->m_rank) && (!a_xg)) {
    fprintf(stderr,"Error in MPIGatherArraynD(): global array is not allocated on root processor.\n");
    return(1);
  }

  /* calculate total size of local domain (w/o a_ghosts) */
  size = 1;
  for (d = 0; d < a_ndims; d++) size *= a_dim_local[d];

  /* create and copy data to send to root process */
  double *buffer = (double*) calloc (size*a_nvars, sizeof(double));
  IERR ArrayCopynD(a_ndims,a_x,buffer,a_dim_local,a_ghosts,0,index,a_nvars); CHECKERR(ierr);

  if (!mpi->m_rank) {
    int proc;
    for (proc = 0; proc < mpi->m_nproc; proc++) {
      int d,done,size;
      /* Find out the domain limits for each process */
      IERR MPILocalDomainLimits(a_ndims,proc,mpi,a_dim_global,is,ie); CHECKERR(ierr);
      size = 1;
      for (d=0; d<a_ndims; d++) {
        size *= (ie[d]-is[d]);
        bounds[d] = ie[d] - is[d];
      }
      if (proc) {
#ifndef serial
        MPI_Status status;
        double *recvbuf = (double*) calloc (size*a_nvars, sizeof(double));
        MPI_Recv(recvbuf,size*a_nvars,MPI_DOUBLE,proc,1902,mpi->m_world,&status);
        int done = 0; _ArraySetValue_(index,a_ndims,0);
        while (!done) {
          int p1; _ArrayIndex1D_(a_ndims,bounds,index,0,p1);
          int p2; _ArrayIndex1DWO_(a_ndims,a_dim_global,index,is,0,p2);
          _ArrayCopy1D_((recvbuf+a_nvars*p1),(a_xg+a_nvars*p2),a_nvars);
          _ArrayIncrementIndex_(a_ndims,bounds,index,done);
        }
        free(recvbuf);
#endif
      } else {
        done = 0; _ArraySetValue_(index,a_ndims,0);
        while (!done) {
          int p1; _ArrayIndex1D_(a_ndims,bounds,index,0,p1);
          int p2; _ArrayIndex1DWO_(a_ndims,a_dim_global,index,is,0,p2);
          _ArrayCopy1D_((buffer+a_nvars*p1),(a_xg+a_nvars*p2),a_nvars);
          _ArrayIncrementIndex_(a_ndims,bounds,index,done);
        }
      }
    }

  } else {
#ifndef serial
    /* Meanwhile, on other processes */
    MPI_Send(buffer,size*a_nvars,MPI_DOUBLE,0,1902,mpi->m_world);
#endif
  }

  free(buffer);
  return(0);
}

/*!
  Gathers the contents of an n-dimensional array, partitioned amongst the MPI ranks, in to a global
  array on rank 0. See documentation of MPIExchangeBoundariesnD() for how the n-dimensional array is
  stored in the memory as a single-index array. This is same as MPIGatherArraynD() but where the
  global array *has* ghost points.

  Notes:
  + The global array must have ghost points (the same number as the local arrays).
  + The global array must be allocated only on rank 0. On other ranks, it must be NULL.

*/
int MPIGatherArraynDwGhosts(
                            int    a_ndims,       /*!< Number of spatial dimensions */
                            void   *a_m,          /*!< MPI object of type #MPIVariables */
                            double *a_xg,         /*!< Global array (preallocated) with ghost points */
                            double *a_x,          /*!< Local array */
                            int    *a_dim_global, /*!< Integer array with elements as global size along each spatial dimension */
                            int    *a_dim_local,  /*!< Integer array with elements as local size along each spatial dimension */
                            int    a_ghosts,      /*!< Number of ghost points */
                            int    a_nvars        /*!< Number of variables (vector components) */
                           )
{
  MPIVariables *mpi = (MPIVariables*) a_m;
  int          d, size;
  _DECLARE_IERR_;

  /* do an exchange to make sure interior/periodic ghost points are filled with consistent values */
  IERR MPIExchangeBoundariesnD(a_ndims, a_nvars, a_dim_local, a_ghosts, mpi, a_x);

  int is[a_ndims], ie[a_ndims], index[a_ndims], bounds[a_ndims];
  int dim_global_wghosts[a_ndims];
  for (d = 0; d < a_ndims; d++) dim_global_wghosts[d] = a_dim_global[d] + 2*a_ghosts;

  /* a_xg should be non-null only on root */
  if (mpi->m_rank && a_xg) {
    fprintf(stderr,"Error in MPIGatherArraynD(): global array exists on non-root processors (rank %d).\n",
            mpi->m_rank);
    return(1);
  }
  if ((!mpi->m_rank) && (!a_xg)) {
    fprintf(stderr,"Error in MPIGatherArraynD(): global array is not allocated on root processor.\n");
    return(1);
  }

  /* calculate total size of local domain (w/ a_ghosts) */
  size = 1;
  for (d = 0; d < a_ndims; d++) size *= (a_dim_local[d]+2*a_ghosts);

  /* create and copy data (incl. ghost points) to send to root process */
  double *buffer = (double*) calloc (size*a_nvars, sizeof(double));
  _ArrayCopy1D_(a_x, buffer, size*a_nvars);

  if (!mpi->m_rank) {
    int proc;
    for (proc = 0; proc < mpi->m_nproc; proc++) {
      int d,done;
      /* Find out the domain limits for each process */
      IERR MPILocalDomainLimits(a_ndims,proc,mpi,a_dim_global,is,ie); CHECKERR(ierr);
      for (d=0; d<a_ndims; d++) {
        bounds[d] = ie[d] - is[d];
      }
      if (proc) {
#ifndef serial
        int size = 1;
        for (d=0; d<a_ndims; d++) {
          size *= (ie[d]-is[d]+2*a_ghosts);
          bounds[d] = ie[d] - is[d];
        }
        MPI_Status status;
        double *recvbuf = (double*) calloc (size*a_nvars, sizeof(double));
        MPI_Recv(recvbuf,size*a_nvars,MPI_DOUBLE,proc,1902,mpi->m_world,&status);
        int done = 0; _ArraySetValue_(index,a_ndims,0);
        while (!done) {
          int p1; _ArrayIndex1D_(a_ndims,bounds,index,a_ghosts,p1);
          int p2; _ArrayIndex1DWO_(a_ndims,a_dim_global,index,is,a_ghosts,p2);
          _ArrayCopy1D_((recvbuf+a_nvars*p1),(a_xg+a_nvars*p2),a_nvars);
          _ArrayIncrementIndex_(a_ndims,bounds,index,done);
        }
        free(recvbuf);
#endif
      } else {
        done = 0; _ArraySetValue_(index,a_ndims,0);
        while (!done) {
          int p1; _ArrayIndex1D_(a_ndims,bounds,index,a_ghosts,p1);
          int p2; _ArrayIndex1DWO_(a_ndims,a_dim_global,index,is,a_ghosts,p2);
          _ArrayCopy1D_((buffer+a_nvars*p1),(a_xg+a_nvars*p2),a_nvars);
          _ArrayIncrementIndex_(a_ndims,bounds,index,done);
        }
      }
    }

  } else {
#ifndef serial
    /* Meanwhile, on other processes */
    MPI_Send(buffer,size*a_nvars,MPI_DOUBLE,0,1902,mpi->m_world);
#endif
  }

  free(buffer);
  return(0);
}
