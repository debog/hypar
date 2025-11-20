/*! @file MPIPartitionArraynD.c
    @brief Partition a global n-dimensional array
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

/*!
  Partitions the contents of a global n-dimensional array on rank 0 (root) to local n-dimensional arrays
  on all the MPI ranks. See documentation of MPIExchangeBoundariesnD() for how the n-dimensional array is
  stored in the memory as a single-index array.

  Notes:
  + The global array must have no ghost points.
  + The global array must be allocated only on rank 0. On other ranks, it must be NULL.
*/
int MPIPartitionArraynD(
                          int     a_ndims,        /*!< Number of spatial dimensions */
                          void    *a_m,           /*!< MPI object of type #MPIVariables */
                          double  *a_xg,          /*!< Global array (preallocated) without ghost points */
                          double  *a_x,           /*!< Local array */
                          int     *a_dim_global,  /*!< Integer array with elements as global size along each spatial dimension */
                          int     *a_dim_local,   /*!< Integer array with elements as local size along each spatial dimension */
                          int     a_ghosts,       /*!< Number of ghost points */
                          int a_nvars /*!< Number of variables (vector components) */
                       )
{
  MPIVariables *mpi = (MPIVariables*) a_m;
  _DECLARE_IERR_;

  int is[a_ndims], ie[a_ndims], index[a_ndims], bounds[a_ndims];

  /* a_xg should be non-null only on root */
  if (mpi->m_rank && a_xg) {
    fprintf(stderr,"Error in MPIPartitionArraynD(): global array exists on non-root processors (rank %d).\n",
            mpi->m_rank);
    return(1);
  }
  if ((!mpi->m_rank) && (!a_xg)) {
    fprintf(stderr,"Error in MPIPartitionArraynD(): global array is not allocated on root processor.\n");
    return(1);
  }

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
      double *buffer = (double*) calloc (size*a_nvars, sizeof(double));
      done = 0; _ArraySetValue_(index,a_ndims,0);
      while (!done) {
        int p1; _ArrayIndex1DWO_(a_ndims,a_dim_global,index,is,0,p1);
        int p2; _ArrayIndex1D_(a_ndims,bounds,index,0,p2);
        _ArrayCopy1D_((a_xg+a_nvars*p1),(buffer+a_nvars*p2),a_nvars);
        _ArrayIncrementIndex_(a_ndims,bounds,index,done);
      }
      if (proc) {
#ifndef serial
        MPI_Send(buffer,size*a_nvars,MPI_DOUBLE,proc,1538,mpi->m_world);
#endif
      } else {
        done = 0; _ArraySetValue_(index,a_ndims,0);
        while (!done) {
          int p1; _ArrayIndex1D_(a_ndims,a_dim_local,index,a_ghosts,p1);
          int p2; _ArrayIndex1D_(a_ndims,a_dim_local,index,0,p2);
          _ArrayCopy1D_((buffer+a_nvars*p2),(a_x+a_nvars*p1),a_nvars);
          _ArrayIncrementIndex_(a_ndims,a_dim_local,index,done);
        }
      }
      free(buffer);
    }

  } else {
#ifndef serial
    /* Meanwhile, on other processes */
    MPI_Status  status;
    int d, done, size;
    size = 1; for (d=0; d<a_ndims; d++) size *= a_dim_local[d];
    double *buffer = (double*) calloc (size*a_nvars, sizeof(double));
    MPI_Recv(buffer,size*a_nvars,MPI_DOUBLE,0,1538,mpi->m_world,&status);
    done = 0; _ArraySetValue_(index,a_ndims,0);
    while (!done) {
      int p1; _ArrayIndex1D_(a_ndims,a_dim_local,index,a_ghosts,p1);
      int p2; _ArrayIndex1D_(a_ndims,a_dim_local,index,0,p2);
      _ArrayCopy1D_((buffer+a_nvars*p2),(a_x+a_nvars*p1),a_nvars);
      _ArrayIncrementIndex_(a_ndims,a_dim_local,index,done);
    }
    free(buffer);
#endif
  }
  return(0);
}

/*!
  Partitions the contents of a global n-dimensional array on rank 0 (root) to local n-dimensional arrays
  on all the MPI ranks. See documentation of MPIExchangeBoundariesnD() for how the n-dimensional array is
  stored in the memory as a single-index array. This is same as MPIPartitionArraynD() but where the global
  array *has* ghost points.

  Notes:
  + The global array must have ghost points (the same number as the local arrays).
  + The global array must be allocated only on rank 0. On other ranks, it must be NULL.
*/
int MPIPartitionArraynDwGhosts(
                                 int     a_ndims,        /*!< Number of spatial dimensions */
                                 void    *a_m,           /*!< MPI object of type #MPIVariables */
                                 double  *a_xg,          /*!< Global array (preallocated) without ghost points */
                                 double  *a_x,           /*!< Local array */
                                 int     *a_dim_global,  /*!< Integer array with elements as global size along each spatial dimension */
                                 int     *a_dim_local,   /*!< Integer array with elements as local size along each spatial dimension */
                                 int     a_ghosts,       /*!< Number of ghost points */
                                 int a_nvars /*!< Number of variables (vector components) */
                              )
{
  MPIVariables *mpi = (MPIVariables*) a_m;
  int d;
  _DECLARE_IERR_;

  int is[a_ndims], ie[a_ndims], index[a_ndims], bounds[a_ndims];
  int dim_global_wghosts[a_ndims];
  for (d = 0; d < a_ndims; d++) dim_global_wghosts[d] = a_dim_global[d] + 2*a_ghosts;

  /* a_xg should be non-null only on root */
  if (mpi->m_rank && a_xg) {
    fprintf(stderr,"Error in MPIPartitionArraynD(): global array exists on non-root processors (rank %d).\n",
            mpi->m_rank);
    return(1);
  }
  if ((!mpi->m_rank) && (!a_xg)) {
    fprintf(stderr,"Error in MPIPartitionArraynD(): global array is not allocated on root processor.\n");
    return(1);
  }

  if (!mpi->m_rank) {
    int proc;
    for (proc = 0; proc < mpi->m_nproc; proc++) {
      int done,size;
      /* Find out the domain limits for each process */
      IERR MPILocalDomainLimits(a_ndims,proc,mpi,a_dim_global,is,ie); CHECKERR(ierr);
      size = 1;
      for (d=0; d<a_ndims; d++) {
        size *= (ie[d]-is[d]+2*a_ghosts);
        bounds[d] = ie[d] - is[d] + 2*a_ghosts;
      }
      double *buffer = (double*) calloc (size*a_nvars, sizeof(double));
      done = 0; _ArraySetValue_(index,a_ndims,0);
      while (!done) {
        int p1; _ArrayIndex1DWO_(a_ndims,dim_global_wghosts,index,is,0,p1);
        int p2; _ArrayIndex1D_(a_ndims,bounds,index,0,p2);
        _ArrayCopy1D_((a_xg+a_nvars*p1),(buffer+a_nvars*p2),a_nvars);
        _ArrayIncrementIndex_(a_ndims,bounds,index,done);
      }
      if (proc) {
#ifndef serial
        MPI_Send(buffer,size*a_nvars,MPI_DOUBLE,proc,1538,mpi->m_world);
#endif
      } else {
        done = 0; _ArraySetValue_(index,a_ndims,0);
        _ArrayCopy1D_(buffer, a_x, size*a_nvars);
      }
      free(buffer);
    }

  } else {
#ifndef serial
    /* Meanwhile, on other processes */
    MPI_Status  status;
    int done, size;
    size = 1; for (d=0; d<a_ndims; d++) size *= (a_dim_local[d]+2*a_ghosts);
    double *buffer = (double*) calloc (size*a_nvars, sizeof(double));
    MPI_Recv(buffer,size*a_nvars,MPI_DOUBLE,0,1538,mpi->m_world,&status);
    _ArrayCopy1D_(buffer, a_x, size*a_nvars);
    free(buffer);
#endif
  }
  return(0);
}
