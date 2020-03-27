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
                          int     ndims,        /*!< Number of spatial dimensions */
                          void    *m,           /*!< MPI object of type #MPIVariables */
                          double  *xg,          /*!< Global array (preallocated) without ghost points */
                          double  *x,           /*!< Local array */
                          int     *dim_global,  /*!< Integer array with elements as global size along each spatial dimension */
                          int     *dim_local,   /*!< Integer array with elements as local size along each spatial dimension */
                          int     ghosts,       /*!< Number of ghost points */
                          int     nvars         /*!< Number of variables (vector components) */
                       )
{
  MPIVariables *mpi = (MPIVariables*) m;
  _DECLARE_IERR_;

  int is[ndims], ie[ndims], index[ndims], bounds[ndims];

  /* xg should be non-null only on root */
  if (mpi->rank && xg) {
    fprintf(stderr,"Error in MPIPartitionArraynD(): global array exists on non-root processors (rank %d).\n",
            mpi->rank);
    return(1);
  }
  if ((!mpi->rank) && (!xg)) {
    fprintf(stderr,"Error in MPIPartitionArraynD(): global array is not allocated on root processor.\n");
    return(1);
  }

  if (!mpi->rank) {
    int proc;
    for (proc = 0; proc < mpi->nproc; proc++) {
      int d,done,size;
      /* Find out the domain limits for each process */
      IERR MPILocalDomainLimits(ndims,proc,mpi,dim_global,is,ie); CHECKERR(ierr);
      size = 1;
      for (d=0; d<ndims; d++) {
        size *= (ie[d]-is[d]);
        bounds[d] = ie[d] - is[d];
      }
      double *buffer = (double*) calloc (size*nvars, sizeof(double));
      done = 0; _ArraySetValue_(index,ndims,0);
      while (!done) {
        int p1; _ArrayIndex1DWO_(ndims,dim_global,index,is,0,p1);
        int p2; _ArrayIndex1D_(ndims,bounds,index,0,p2);
        _ArrayCopy1D_((xg+nvars*p1),(buffer+nvars*p2),nvars);
        _ArrayIncrementIndex_(ndims,bounds,index,done);
      }
      if (proc) {
#ifndef serial
        MPI_Send(buffer,size*nvars,MPI_DOUBLE,proc,1538,mpi->world);
#endif
      } else {
        done = 0; _ArraySetValue_(index,ndims,0);
        while (!done) {
          int p1; _ArrayIndex1D_(ndims,dim_local,index,ghosts,p1);
          int p2; _ArrayIndex1D_(ndims,dim_local,index,0,p2);
          _ArrayCopy1D_((buffer+nvars*p2),(x+nvars*p1),nvars);
          _ArrayIncrementIndex_(ndims,dim_local,index,done);
        }
      }
      free(buffer);
    }

  } else {
#ifndef serial
    /* Meanwhile, on other processes */
    MPI_Status  status;
    int d, done, size;
    size = 1; for (d=0; d<ndims; d++) size *= dim_local[d];
    double *buffer = (double*) calloc (size*nvars, sizeof(double));
    MPI_Recv(buffer,size*nvars,MPI_DOUBLE,0,1538,mpi->world,&status);
    done = 0; _ArraySetValue_(index,ndims,0);
    while (!done) {
      int p1; _ArrayIndex1D_(ndims,dim_local,index,ghosts,p1);
      int p2; _ArrayIndex1D_(ndims,dim_local,index,0,p2);
      _ArrayCopy1D_((buffer+nvars*p2),(x+nvars*p1),nvars);
      _ArrayIncrementIndex_(ndims,dim_local,index,done);
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
                                 int     ndims,        /*!< Number of spatial dimensions */
                                 void    *m,           /*!< MPI object of type #MPIVariables */
                                 double  *xg,          /*!< Global array (preallocated) without ghost points */
                                 double  *x,           /*!< Local array */
                                 int     *dim_global,  /*!< Integer array with elements as global size along each spatial dimension */
                                 int     *dim_local,   /*!< Integer array with elements as local size along each spatial dimension */
                                 int     ghosts,       /*!< Number of ghost points */
                                 int     nvars         /*!< Number of variables (vector components) */
                              )
{
  MPIVariables *mpi = (MPIVariables*) m;
  int d;
  _DECLARE_IERR_;

  int is[ndims], ie[ndims], index[ndims], bounds[ndims];
  int dim_global_wghosts[ndims];
  for (d = 0; d < ndims; d++) dim_global_wghosts[d] = dim_global[d] + 2*ghosts;

  /* xg should be non-null only on root */
  if (mpi->rank && xg) {
    fprintf(stderr,"Error in MPIPartitionArraynD(): global array exists on non-root processors (rank %d).\n",
            mpi->rank);
    return(1);
  }
  if ((!mpi->rank) && (!xg)) {
    fprintf(stderr,"Error in MPIPartitionArraynD(): global array is not allocated on root processor.\n");
    return(1);
  }

  if (!mpi->rank) {
    int proc;
    for (proc = 0; proc < mpi->nproc; proc++) {
      int done,size;
      /* Find out the domain limits for each process */
      IERR MPILocalDomainLimits(ndims,proc,mpi,dim_global,is,ie); CHECKERR(ierr);
      size = 1;
      for (d=0; d<ndims; d++) {
        size *= (ie[d]-is[d]+2*ghosts);
        bounds[d] = ie[d] - is[d] + 2*ghosts;
      }
      double *buffer = (double*) calloc (size*nvars, sizeof(double));
      done = 0; _ArraySetValue_(index,ndims,0);
      while (!done) {
        int p1; _ArrayIndex1DWO_(ndims,dim_global_wghosts,index,is,0,p1);
        int p2; _ArrayIndex1D_(ndims,bounds,index,0,p2);
        _ArrayCopy1D_((xg+nvars*p1),(buffer+nvars*p2),nvars);
        _ArrayIncrementIndex_(ndims,bounds,index,done);
      }
      if (proc) {
#ifndef serial
        MPI_Send(buffer,size*nvars,MPI_DOUBLE,proc,1538,mpi->world);
#endif
      } else {
        done = 0; _ArraySetValue_(index,ndims,0);
        _ArrayCopy1D_(buffer, x, size*nvars);
      }
      free(buffer);
    }

  } else {
#ifndef serial
    /* Meanwhile, on other processes */
    MPI_Status  status;
    int done, size;
    size = 1; for (d=0; d<ndims; d++) size *= (dim_local[d]+2*ghosts);
    double *buffer = (double*) calloc (size*nvars, sizeof(double));
    MPI_Recv(buffer,size*nvars,MPI_DOUBLE,0,1538,mpi->world,&status);
    _ArrayCopy1D_(buffer, x, size*nvars);
    free(buffer);
#endif
  }
  return(0);
}
